import os
import sys
import time
import re
import csv
from datetime import datetime
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def run_step1b_final_dual_output():
    start_time = time.time()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1A")
    output_path = os.path.join(project_root, "Step_Outputs", "Phase1B")
    raw_out = os.path.join(output_path, "Raw_Stability")
    filt_out = os.path.join(output_path, "Filtered_Stability")
    os.makedirs(raw_out, exist_ok=True)
    os.makedirs(filt_out, exist_ok=True)
    
    if not os.path.exists(input_folder):
        print("\n[ERROR] Phase 1A directory not found.")
        return

    fasta_files = [f for f in os.listdir(input_folder) if f.endswith(".fasta")]
    total_files = len(fasta_files)

    print("\n" + "="*80)
    print(f"{'PHASE 1B: STABILITY ANALYSIS (PROTPARAM)':^80}")
    print("="*80)

    raw_list, filtered_list = [], []
    sys.stdout.write("[PROCESS] Initiating ExPASy analytical engine...\n\n")
    
    for i, filename in enumerate(sorted(fasta_files)):
        file_path = os.path.join(input_folder, filename)
        variant_id = filename.replace(".fasta", "")
        target_group = filename.split('_Var')[0]
        
        elapsed = format_time(time.time() - start_time)
        sys.stdout.write(f"\r[ EVAL ] Record {i+1:03d}/{total_files:03d} | Target: {target_group:<12} | Elapsed: {elapsed}")
        sys.stdout.flush()
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', str(record.seq).upper())
                if not clean_seq: continue

                analysis = ProteinAnalysis(clean_seq)
                mw = analysis.molecular_weight() / 1000
                gravy = analysis.gravy()
                idx = analysis.instability_index()
                
                status = "STABLE" if idx < 40 else "UNSTABLE"
                data_row = {"Variant_ID": variant_id, "Target": target_group, "MW_kDa": round(mw, 2), "GRAVY": round(gravy, 2), "Instability_Index": round(idx, 2), "Status": status}

                raw_list.append(data_row)
                if idx < 40:
                    filtered_list.append(data_row)
                    # [PATCH]: Save the surviving FASTA for Phase 1C
                    with open(os.path.join(filt_out, filename), "w") as out_fasta:
                        SeqIO.write(record, out_fasta, "fasta")
        except Exception as e:
            continue

    print() 
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M")
    raw_csv = os.path.join(raw_out, f"Phase1B_Raw_Full_{timestamp}.csv")
    filt_csv = os.path.join(filt_out, f"Phase1B_Filtered_Stable_{timestamp}.csv")
    FIELD_NAMES = ["Variant_ID", "Target", "MW_kDa", "GRAVY", "Instability_Index", "Status"]

    if raw_list:
        with open(raw_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=FIELD_NAMES)
            writer.writeheader()
            writer.writerows(raw_list)
        if filtered_list:
            with open(filt_csv, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=FIELD_NAMES)
                writer.writeheader()
                writer.writerows(filtered_list)

    print("\n" + "="*80)
    print(f"[SUCCESS] Stable Candidates Passed : {len(filtered_list)} / {len(raw_list)}")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_step1b_final_dual_output()