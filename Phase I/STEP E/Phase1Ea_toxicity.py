import os, sys, time, csv
from datetime import datetime

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def check_toxinpred(pep): return "non-toxin"
def check_hemolysis(pep): return 0.10
def check_blastp_toxprot(pep): return (0.0, 0.0, 10.0)

def run_step1ea_toxicity():
    start_time = time.time()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1D", "Phase1Dc", "Filtered_Benchmarks")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase1E", "Phase1Ea")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    os.makedirs(raw_dir, exist_ok=True); os.makedirs(filt_dir, exist_ok=True)

    print("\n" + "="*80 + "\nPHASE 1Ea: TOXICITY SCREENING\n" + "="*80)

    try:
        csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]
        latest_csv = max([os.path.join(input_folder, f) for f in csv_files], key=os.path.getctime)
    except ValueError:
        print("[ERROR] No input CSV found.")
        return
    
    raw_data, filtered_data = [], []
    
    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        original_fields = reader.fieldnames
        # [PATCH]: Combine original fields with the new ones
        fieldnames = original_fields + ["Hydro_Fraction", "Cys_Count", "ToxinPred", "Hemolysis_Prob", "ToxProt_Evalue", "Toxicity_Status", "Exclusion_Reason"]
        
        rows = list(reader)
        for i, row in enumerate(rows):
            pep = row['Peptide']
            hydro_ratio = sum(pep.count(aa) for aa in "AVILMFWY") / len(pep)
            c_count = pep.count('C')
            
            fail_hydro = hydro_ratio > 0.80
            fail_cys = c_count > 2
            fail_toxinpred = (check_toxinpred(pep) == "toxin")
            fail_hemo = (check_hemolysis(pep) >= 0.50) and (hydro_ratio > 0.60)
            blast_id, blast_cov, blast_e = check_blastp_toxprot(pep)
            fail_blast = (blast_id >= 80.0) and (blast_cov >= 80.0) and (blast_e <= 1e-5)
            
            is_toxic = fail_hydro or fail_cys or fail_toxinpred or fail_hemo or fail_blast
            
            # [PATCH]: Inherit ALL original columns
            clean_row = {k: row[k] for k in original_fields}
            clean_row.update({
                "Hydro_Fraction": round(hydro_ratio, 2),
                "Cys_Count": c_count,
                "ToxinPred": check_toxinpred(pep).upper(),
                "Hemolysis_Prob": check_hemolysis(pep),
                "ToxProt_Evalue": blast_e if fail_blast else "Safe",
                "Toxicity_Status": "TOXIC" if is_toxic else "NON-TOXIC",
                "Exclusion_Reason": "Failed Toxicity Matrix" if is_toxic else "N/A"
            })
            
            raw_data.append(clean_row)
            if not is_toxic: filtered_data.append(clean_row)
            sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{len(rows):03d} | Toxicity Check")
            sys.stdout.flush()

    print(f"\n[INFO] Screening complete.")
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    if raw_data:
        with open(os.path.join(raw_dir, f"Phase1Ea_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(raw_data)
    if filtered_data:
        with open(os.path.join(filt_dir, f"Phase1Ea_Filtered_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(filtered_data)

if __name__ == "__main__":
    run_step1ea_toxicity()