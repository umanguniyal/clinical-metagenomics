import argparse
import json
import os

def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        try:
            return json.load(f)
        except:
            return None

def extract_entities(report):
    species = set()
    genes = set()
    if not report:
        return species, genes
    
    # Extract species
    for tax_key in ["kraken2", "centrifuge"]:
        if tax_key in report and "taxonomy" in report[tax_key]:
            for tax in report[tax_key]["taxonomy"]:
                if "name" in tax:
                    species.add(tax["name"])

    # Extract AMR genes
    amr_section = report.get("amr", {})
    for tool in ["abricate", "amrfinder", "kma"]:
        if tool in amr_section:
            for item in amr_section[tool]:
                if "gene" in item:
                    genes.add(item["gene"])
                elif "Gene" in item:
                    genes.add(item["Gene"])
                    
    return species, genes

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient_json", required=True)
    ap.add_argument("--control_json", required=True)
    ap.add_argument("--out_json", required=True)
    args = ap.parse_args()

    patient = load_json(args.patient_json)
    control = load_json(args.control_json)

    if not patient:
        print(f"Patient JSON not found at {args.patient_json}")
        return

    if not control:
        with open(args.out_json, "w") as f:
            json.dump(patient, f, indent=2)
        return

    pat_species, pat_genes = extract_entities(patient)
    ctrl_species, ctrl_genes = extract_entities(control)

    shared_species = pat_species.intersection(ctrl_species)
    shared_genes = pat_genes.intersection(ctrl_genes)

    flags = []
    if shared_species:
        flags.append(f"Organisms detected in Negative Control: {', '.join(sorted(shared_species))}")
    if shared_genes:
        flags.append(f"AMR Genes detected in Negative Control: {', '.join(sorted(shared_genes))}")

    if flags:
        patient["contamination_flags"] = flags
    else:
        patient["contamination_flags"] = []

    with open(args.out_json, "w") as f:
        json.dump(patient, f, indent=2)

if __name__ == "__main__":
    main()
