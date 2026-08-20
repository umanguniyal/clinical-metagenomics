import os
import glob
import re

RULES_DIR = r"d:\WSL\Ubuntu\rootfs\root\UMANG\SRR\clinical\pipeline_smk\rules"
smk_files = glob.glob(os.path.join(RULES_DIR, "*.smk"))

for smk_file in smk_files:
    with open(smk_file, "r", encoding="utf-8") as f:
        content = f.read()

    new_content = []
    lines = content.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Detect rule definition
        match_rule = re.match(r'^\s*rule\s+([A-Za-z0-9_]+):', line)
        if match_rule:
            rule_name = match_rule.group(1)
            new_content.append(line)
            
            # Find the extent of this rule
            j = i + 1
            rule_lines = []
            while j < len(lines):
                # Next rule or top level block
                # Watch out for indented lines which are part of the rule
                if re.match(r'^\s*rule\s+', lines[j]) or (re.match(r'^[a-zA-Z_]', lines[j]) and not lines[j].startswith('if ') and not lines[j].startswith('elif ') and not lines[j].startswith('else:')):
                    break
                rule_lines.append(lines[j])
                j += 1
                
            in_log = False
            in_output = False
            log_end_idx = -1
            output_end_idx = -1
            
            pattern_line = ""
            
            for k, rline in enumerate(rule_lines):
                if re.match(r'^\s+output:', rline):
                    in_output = True
                    in_log = False
                elif re.match(r'^\s+log:', rline):
                    in_log = True
                    in_output = False
                elif re.match(r'^\s+[a-zA-Z0-9_]+:', rline):
                    in_log = False
                    in_output = False
                    
                if in_output or in_log:
                    if "RESULTS_DIR" in rline:
                        pattern_line = rline
                        
                if in_log:
                    log_end_idx = k
                if in_output:
                    output_end_idx = k
                    
            insert_idx = log_end_idx if log_end_idx != -1 else output_end_idx
                
            if not pattern_line:
                for rline in rule_lines:
                    if "RESULTS_DIR" in rline:
                        pattern_line = rline
                        break
                        
            platform_str = ""
            sample_str = ""
            if '"{platform}"' in pattern_line:
                platform_str = '"{platform}", '
            elif '"illumina"' in pattern_line:
                platform_str = '"illumina", '
            elif '"nanopore"' in pattern_line:
                platform_str = '"nanopore", '
                
            if '"{sample}"' in pattern_line:
                sample_str = '"{sample}", '

            bench_path = f'        os.path.join(RESULTS_DIR, {platform_str}{sample_str}"benchmarks", "{rule_name}.txt")'
            
            for k, rline in enumerate(rule_lines):
                new_content.append(rline)
                if k == insert_idx:
                    new_content.append("    benchmark:")
                    new_content.append(bench_path)
            
            i = j
        else:
            new_content.append(line)
            i += 1

    with open(smk_file, "w", encoding="utf-8") as f:
        f.write('\n'.join(new_content))
