import argparse
import json
import os
import pandas as pd

from ml_pathogen_amr import extract_features

def evaluate_threat(row, bartlett_species):
    tags = []
    
    # Check Bartlett
    is_known_pathogen = False
    org_name = str(row["organism"])
    
    # Simple check if organism name contains any Bartlett species
    for sp in bartlett_species:
        if isinstance(sp, str) and sp.lower() in org_name.lower():
            is_known_pathogen = True
            break
            
    if is_known_pathogen:
        tags.append("[KNOWN PATHOGEN]")
        
    # Check Abundance
    if row.get("rel_abundance", 0) >= 5.0:
        tags.append("[HIGH ABUNDANCE]")
        
    # Check Coverage/Quality
    if row.get("completeness", 0) >= 70.0 and row.get("contamination", 100) <= 10.0:
        tags.append("[HIGH-QUALITY BIN]")
        
    # Check AMR Precision
    has_binned_amr = row.get("n_amr_genes", 0) >= 1
    has_kma_amr = row.get("kma_n_hits", 0) >= 1
    max_id = row.get("max_amr_identity", 0)
    
    if max_id >= 90.0 and has_binned_amr:
        tags.append("[CONFIRMED BINNED AMR]")
    elif has_kma_amr and not has_binned_amr:
        tags.append("[STATISTICAL KMA AMR]")
        
    # Assign Final Threat Tier
    has_known = "[KNOWN PATHOGEN]" in tags
    has_conf_amr = "[CONFIRMED BINNED AMR]" in tags
    has_hq_bin = "[HIGH-QUALITY BIN]" in tags
    has_high_abund = "[HIGH ABUNDANCE]" in tags
    
    if has_known and has_conf_amr and has_hq_bin:
        tier = "CRITICAL"
    elif has_known and (has_conf_amr or has_high_abund):
        tier = "HIGH"
    elif (has_conf_amr and not has_known) or (has_known and not has_conf_amr and not has_kma_amr):
        tier = "MODERATE"
    else:
        tier = "LOW"
        
    return tier, tags

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report_json", required=True)
    ap.add_argument("--bartlett_xlsx", required=True)
    ap.add_argument("--viral_xlsx", default="")
    ap.add_argument("--fungal_xlsx", default="")
    ap.add_argument("--out_txt", required=True)
    args = ap.parse_args()
    
    COMMON_SYNONYMS = {
        "Clostridium difficile": "Clostridioides difficile",
        "Propionibacterium acnes": "Cutibacterium acnes",
        "Campylobacter pylori": "Helicobacter pylori",
        "Pneumocystis carinii": "Pneumocystis jirovecii",
    }
    
    with open(args.report_json, "r") as f:
        report = json.load(f)
        
    # Load Bartlett DB
    try:
        bartlett_df = pd.read_excel(args.bartlett_xlsx)
        
        # Combine Genus and Species if both exist
        if "Genus" in bartlett_df.columns and "Species" in bartlett_df.columns:
            bartlett_species = (bartlett_df["Genus"].astype(str).str.strip() + " " + 
                              bartlett_df["Species"].astype(str).str.strip()).dropna().unique().tolist()
        else:
            col = "Species" if "Species" in bartlett_df.columns else bartlett_df.columns[0]
            bartlett_species = bartlett_df[col].dropna().unique().tolist()
            
        # If there's a synonym column in Bartlett, load it
        for col in ["Previous name", "Synonym", "Synonyms"]:
            if col in bartlett_df.columns:
                bartlett_species.extend(bartlett_df[col].dropna().astype(str).unique().tolist())
    except Exception as e:
        print(f"Warning: Could not read Bartlett database properly: {e}")
        bartlett_species = []

    # Load Viral pathogen DB
    if args.viral_xlsx and os.path.exists(args.viral_xlsx):
        try:
            vdf = pd.read_excel(args.viral_xlsx, engine="openpyxl")
            for col in ["ICTV species", "Common name", "Synonyms"]:
                if col in vdf.columns:
                    bartlett_species.extend(vdf[col].dropna().astype(str).unique().tolist())
            print(f"Loaded {len(vdf)} viral pathogen entries.")
        except Exception as e:
            print(f"Warning: Could not read viral pathogen database: {e}")

    # Load Fungal pathogen DB (CSV disguised as .xlsx, filter by human.pathogen)
    if args.fungal_xlsx and os.path.exists(args.fungal_xlsx):
        try:
            fdf = pd.read_csv(args.fungal_xlsx)
            if "human.pathogen" in fdf.columns:
                fdf = fdf[fdf["human.pathogen"] == True]
            if "Species" in fdf.columns:
                bartlett_species.extend(fdf["Species"].dropna().astype(str).unique().tolist())
            if "atlas.synonym.name" in fdf.columns:
                bartlett_species.extend(fdf["atlas.synonym.name"].dropna().astype(str).unique().tolist())
            print(f"Loaded {len(fdf)} human-pathogenic fungal entries.")
        except Exception as e:
            print(f"Warning: Could not read fungal pathogen database: {e}")

    # Add hardcoded synonyms
    for old_name, new_name in COMMON_SYNONYMS.items():
        if new_name in bartlett_species:
            bartlett_species.append(old_name)

    # Deduplicate
    bartlett_species = list(set([s.strip() for s in bartlett_species if str(s).strip() != ""]))
        
    # Parse organisms
    df = extract_features(report)
    
    list1_susceptible = []
    list2_amr_threats = []
    
    if not df.empty:
        for _, row in df.iterrows():
            tier, tags = evaluate_threat(row, bartlett_species)
            
            org_data = {
                "organism": row["organism"],
                "rel_abundance": row["rel_abundance"],
                "tier": tier,
                "tags": tags,
                "amr_genes": row.get("amr_gene_names", "")
            }
            
            n_amr = row.get("n_amr_genes", 0)
            
            # List 1: Susceptible Pathogens
            if (row.get("rel_abundance", 0) > 1.0 or row.get("completeness", 0) > 50) and n_amr == 0:
                list1_susceptible.append(org_data)
                
            # List 2: AMR-Associated Threats
            if n_amr >= 1 and (row.get("rel_abundance", 0) > 0.1 or row.get("completeness", 0) > 20):
                list2_amr_threats.append(org_data)
                
    # Sort lists by: 1) Known Pathogen status, 2) tier criticality, 3) relative abundance
    tier_order = {"CRITICAL": 0, "HIGH": 1, "MODERATE": 2, "LOW": 3}
    
    def triage_sort_key(x):
        is_known = 0 if "[KNOWN PATHOGEN]" in x["tags"] else 1
        tier = tier_order.get(x["tier"], 4)
        abund = -x["rel_abundance"]
        return (is_known, tier, abund)
        
    list1_susceptible.sort(key=triage_sort_key)
    list2_amr_threats.sort(key=triage_sort_key)
    
    # Inject into JSON
    report["bioinfo_triage_lists"] = {
        "susceptible_pathogens": list1_susceptible,
        "amr_associated_threats": list2_amr_threats
    }
    
    with open(args.report_json, "w") as f:
        json.dump(report, f, indent=2)
        
    # Append to TXT report
    with open(args.out_txt, "a") as f:
        f.write("\n\n=== Clinical Threat Triage (Bioinformatician) ===\n")
        
        f.write("\n--- Susceptible Pathogen Candidates ---\n")
        if not list1_susceptible:
            f.write("None detected.\n")
        for item in list1_susceptible:
            tags_str = " ".join(item["tags"])
            f.write(f"[{item['tier']}] {item['organism']} (Abundance: {item['rel_abundance']:.1f}%)\n")
            if tags_str:
                f.write(f"    Tags: {tags_str}\n")
                
        f.write("\n--- AMR-Associated Threats ---\n")
        if not list2_amr_threats:
            f.write("None detected.\n")
        for item in list2_amr_threats:
            tags_str = " ".join(item["tags"])
            f.write(f"[{item['tier']}] {item['organism']} (Abundance: {item['rel_abundance']:.1f}%)\n")
            if tags_str:
                f.write(f"    Tags: {tags_str}\n")
            if item["amr_genes"]:
                f.write(f"    AMR Genes: {item['amr_genes']}\n")

if __name__ == "__main__":
    main()
