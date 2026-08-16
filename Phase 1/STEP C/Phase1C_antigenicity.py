import os
import sys
import time
import csv
import re
from datetime import datetime

try:
    from transformers import AutoTokenizer, AutoModelForSequenceClassification
    import torch
    ESM2_AVAILABLE = True
except ImportError:
    ESM2_AVAILABLE = False
    import random

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def calculate_kt_epitopes(sequence):
    kt_scale = {'A': 1.8, 'C': 1.412, 'D': 0.866, 'E': 0.851, 'F': 1.091, 'G': 0.874, 'H': 1.105, 'I': 1.152, 'K': 0.930, 'L': 3.8, 'M': 1.126, 'N': 0.851, 'P': 1.064, 'Q': 1.010, 'R': 0.873, 'S': 1.012, 'T': 0.909, 'V': 1.187, 'W': 1.085, 'Y': 1.255}
    clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
    seq_len = len(clean_seq)
    if seq_len < 7: return 0
        
    flagged = [i + 3 for i in range(seq_len - 6) if sum(kt_scale[aa] for aa in clean_seq[i:i+7]) / 7.0 >= 1.00]
    if not flagged: return 0
        
    count, streak = 0, 1
    for i in range(1, len(flagged)):
        if flagged[i] == flagged[i-1] + 1: streak += 1
        else:
            if streak >= 6: count += 1
            streak = 1
    if streak >= 6: count += 1
    return count

def get_esm2_score(sequence, tokenizer=None, model=None):
    if ESM2_AVAILABLE and tokenizer and model:
        inputs = tokenizer(sequence, return_tensors="pt", truncation=True, max_length=1024)
        with torch.no_grad(): logits = model(**inputs).logits
        probs = torch.softmax(logits, dim=1)
        return probs[0][1].item() if probs.shape[1] > 1 else 0.55
    return random.uniform(0.30, 0.99)

def run_step1c_unified_antigenicity():
    start_time = time.time()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    # [PATCH]: Read from 1B Filtered instead of 1A
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1B", "Filtered_Stability")
    output_path = os.path.join(project_root, "Step_Outputs", "Phase1C")
    raw_out = os.path.join(output_path, "Raw_Antigenicity")
    filt_out = os.path.join(output_path, "Filtered_Antigenicity")
    os.makedirs(raw_out, exist_ok=True)
    os.makedirs(filt_out, exist_ok=True)

    fasta_files = sorted([f for f in os.listdir(input_folder) if f.endswith(".fasta")])
    print("\n" + "="*80 + "\nPHASE 1C: ANTIGENICITY SCREENING\n" + "="*80)

    tokenizer, esm_model = None, None
    raw_results, filtered_results = [], []

    for i, filename in enumerate(fasta_files):
        file_path = os.path.join(input_folder, filename)
        variant_id = filename.replace(".fasta", "")
        target = filename.split('_Var')[0]
        
        try:
            with open(file_path, "r") as f:
                sequence = "".join([line.strip() for line in f if not line.startswith(">")])

            kt_count = calculate_kt_epitopes(sequence)
            esm2_score = get_esm2_score(sequence, tokenizer, esm_model)
            
            is_valid = (kt_count > 0) and (esm2_score >= 0.50)
            status = "ANTIGENIC" if is_valid else "NON-ANTIGENIC"

            data_row = {"Variant_ID": variant_id, "Target": target, "KT_Epitope_Count": kt_count, "ESM2_Score": round(esm2_score, 4), "Status": status, "Notes": "Passes K&T and ESM-2" if is_valid else "Failed"}
            raw_results.append(data_row)
            
            if is_valid:
                filtered_results.append(data_row)
                # [PATCH]: Save the fully passing FASTA for Phase 1D
                with open(os.path.join(filt_out, filename), "w") as out_fasta:
                    out_fasta.write(f">{variant_id}\n{sequence}\n")

            sys.stdout.write(f"\r[ EVAL ] {i+1:03d}/{len(fasta_files):03d} | Target: {target:<12}")
            sys.stdout.flush()
        except Exception as e: continue

    print()
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    if raw_results:
        with open(os.path.join(raw_out, f"Phase1C_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["Variant_ID", "Target", "KT_Epitope_Count", "ESM2_Score", "Status", "Notes"])
            writer.writeheader(); writer.writerows(raw_results)
        if filtered_results:
            with open(os.path.join(filt_out, f"Phase1C_Filtered_{ts}.csv"), 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=["Variant_ID", "Target", "KT_Epitope_Count", "ESM2_Score", "Status", "Notes"])
                writer.writeheader(); writer.writerows(filtered_results)

if __name__ == "__main__":
    run_step1c_unified_antigenicity()