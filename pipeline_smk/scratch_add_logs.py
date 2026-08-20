import os
import glob
import re

def process_smk(file_path):
    with open(file_path, "r") as f:
        content = f.read()

    rules = re.split(r'^(rule\s+\w+:)', content, flags=re.MULTILINE)
    if len(rules) == 1:
        return

    out_parts = [rules[0]]
    for i in range(1, len(rules), 2):
        rule_decl = rules[i]
        rule_name = re.match(r'rule\s+(\w+):', rule_decl).group(1)
        rule_body = rules[i+1]
        
        # Check if log already exists
        if re.search(r'^\s+log:', rule_body, flags=re.MULTILINE):
            out_parts.append(rule_decl + rule_body)
            continue
            
        # Determine if it's a sample-level rule or cohort
        has_sample = "{sample}" in rule_body
        is_illumina = "illumina" in rule_path or "_illumina" in rule_name or "illumina" in rule_body
        is_nanopore = "nanopore" in rule_path or "_nanopore" in rule_name or "nanopore" in rule_body
        
        # Simple heuristic for platform
        if is_illumina and not is_nanopore:
            plat = '"illumina"'
        elif is_nanopore and not is_illumina:
            plat = '"nanopore"'
        else:
            plat = 'PLATFORM'
            
        if has_sample:
            log_path = f'os.path.join(RESULTS_DIR, {plat}, "{{sample}}", "logs", "{rule_name}.log")'
        else:
            log_path = f'os.path.join(RESULTS_DIR, "logs", "{rule_name}.log")'
            
        # Insert log directive before conda or threads or shell
        insert_idx = -1
        for keyword in ["conda:", "threads:", "shell:", "run:"]:
            match = re.search(r'^(\s+)' + keyword, rule_body, flags=re.MULTILINE)
            if match:
                insert_idx = match.start()
                indent = match.group(1)
                break
                
        if insert_idx != -1:
            log_str = f"{indent}log:\n{indent}    {log_path}\n"
            rule_body = rule_body[:insert_idx] + log_str + rule_body[insert_idx:]
            
            # Now try to append > {log} 2>&1 to shell block
            shell_match = re.search(r'^\s+shell:\s*(r?"""|r?\'\'\')(.*?)\1', rule_body, flags=re.MULTILINE|re.DOTALL)
            if shell_match:
                shell_content = shell_match.group(2)
                # If there's already redirection, skip
                if "> {log}" not in shell_content:
                    # Append it right before the end
                    new_shell = shell_content.rstrip() + " > {log} 2>&1\n        "
                    rule_body = rule_body[:shell_match.start(2)] + new_shell + rule_body[shell_match.end(2):]
        
        out_parts.append(rule_decl + rule_body)
        
    with open(file_path, "w") as f:
        f.write("".join(out_parts))

if __name__ == "__main__":
    for rule_path in glob.glob("d:/WSL/Ubuntu/rootfs/root/UMANG/SRR/clinical/pipeline_smk/rules/*.smk"):
        process_smk(rule_path)
