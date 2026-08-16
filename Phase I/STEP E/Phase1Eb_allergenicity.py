import os, sys, time, csv
from datetime import datetime

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def run_step1eb_allergenicity():
    start_time = time.time()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1E", "Phase1Ea", "Filtered")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase1E", "Phase1Eb")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    os.makedirs(raw_dir, exist_ok=True); os.makedirs(filt_dir, exist_ok=True)

    print("\n" + "="*80 + "\nPHASE 1Eb: ALLERGENICITY SCREENING\n" + "="*80)

    files = sorted([f for f in os.listdir(input_folder) if f.endswith(".csv")])
    if not files: return
    latest_csv = os.path.join(input_folder, files[-1])
    
    raw_data, filtered_data = [], []
    
    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        original_fields = reader.fieldnames
        # [PATCH]: Combine original fields with the new ones
        fieldnames = original_fields + ["QN_Ratio", "Surface_Charge", "Allergen_Status"]
        
        rows = list(reader)
        for i, row in enumerate(rows):
            pep = row['Peptide']
            qn_ratio = (pep.count('Q') + pep.count('N')) / len(pep)
            charge_count = sum(pep.count(aa) for aa in "DEHK")
            is_allergen = (qn_ratio > 0.3) or (charge_count > 4)
            
            # [PATCH]: Inherit ALL original columns
            clean_row = {k: row[k] for k in original_fields}
            clean_row.update({
                "QN_Ratio": round(qn_ratio, 2),
                "Surface_Charge": charge_count,
                "Allergen_Status": "ALLERGEN" if is_allergen else "NON-ALLERGEN"
            })
            
            raw_data.append(clean_row)
            if not is_allergen: filtered_data.append(clean_row)
            sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{len(rows):03d} | Allergen Check")
            sys.stdout.flush()

    print(f"\n[INFO] Screening complete.")
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    if raw_data:
        with open(os.path.join(raw_dir, f"Phase1Eb_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(raw_data)
    if filtered_data:
        with open(os.path.join(filt_dir, f"Phase1Eb_Filtered_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(filtered_data)

if __name__ == "__main__":
    run_step1eb_allergenicity()