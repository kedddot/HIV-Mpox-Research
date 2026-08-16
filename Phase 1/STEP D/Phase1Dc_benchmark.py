import os
import csv
import re
import time
from datetime import datetime
# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================
def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"
# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================
def run_step1dc_conservancy_benchmark():
    start_time = time.time()
    
    # DYNAMIC PATHING
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    # Inputs: Use latest from 1Db and original FASTA library
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1D", "Phase1Db")
    fasta_folder = os.path.join(project_root, "Step_Outputs", "Phase1C", "Filtered_Antigenicity")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase1D", "Phase1Dc")
    
    raw_dir = os.path.join(output_base, "Raw_Conservancy")
    filt_dir = os.path.join(output_base, "Filtered_Benchmarks")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)
    print("\n" + "="*80)
    print(f"{'PHASE 1Dc: CONSERVANCY & VARIANT BENCHMARKING':^80}")
    print("="*80)
    # 1. LOAD PROTEIN LIBRARY
    print(f"[PROCESS] Loading viral protein library from {os.path.basename(fasta_folder)}...")
    library = {}
    for f in os.listdir(fasta_folder):
        if f.endswith(".fasta"):
            with open(os.path.join(fasta_folder, f), "r") as file:
                seq = "".join([l.strip() for l in file if not l.startswith(">")])
                library[f] = {"seq": re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper()), "target": f.split('_Var')[0]}
    # 2. LOAD ELITE BINDERS
    db_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]
    if not db_files:
        print("[ERROR] No input CSV found in Phase 1Db.")
        return
    latest_db = os.path.join(input_folder, sorted(db_files)[-1])
    
    print(f"[PROCESS] Analyzing conservancy against {len(library)} variants...")
    final_results = []
    
    with open(latest_db, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pep = row['Peptide']
            target = row['Target']
            relevant_vars = [v['seq'] for v in library.values() if v['target'] == target]
            total_relevant = len(relevant_vars) if relevant_vars else 1
            hits = sum(1 for s in relevant_vars if pep in s)
            
            row['Conservancy'] = round((hits / total_relevant) * 100, 2)
            row['Hit_Ratio'] = f"{hits}/{total_relevant}"
            final_results.append(row)
    # 3. EXPORT
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    raw_path = os.path.join(raw_dir, f"Phase1Dc_Raw_Full_{ts}.csv")
    
    with open(raw_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=final_results[0].keys())
        writer.writeheader()
        writer.writerows(final_results)
    # BENCHMARKS (Thresholds)
    # NOTE: methods specify strict inequalities ("> 20%", "> 50%") for the minimal and
    # highly-acceptable bands, with 100% as the (inclusive, since it's a ceiling) ideal
    # goal. The previous version used >= for all three, which incorrectly included
    # peptides sitting exactly at 20.00% or 50.00% in bands the text excludes them from.
    print("-" * 80)
    print(f"{'BENCHMARK THRESHOLD RESULTS':^80}")
    thresholds = [
        (20, "Minimal (> 20%)", lambda c: c > 20),
        (50, "Highly Acceptable (> 50%)", lambda c: c > 50),
        (100, "Ideal (100%)", lambda c: c >= 100),
    ]
    for t, label, test in thresholds:
        filtered = [r for r in final_results if test(float(r['Conservancy']))]
        if filtered:
            f_path = os.path.join(filt_dir, f"Phase1Dc_Min_{t}pct_{ts}.csv")
            with open(f_path, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=final_results[0].keys())
                writer.writeheader()
                writer.writerows(filtered)
            print(f"[SUCCESS] {label}: {len(filtered)} epitopes found.")
            # Per-Type breakdown -- Phase1Db now emits MHC-I, MHC-II, and B-cell rows
            # together, so a flat count here would otherwise hide the split.
            if filtered and "Type" in filtered[0]:
                type_counts = {}
                for r in filtered:
                    type_counts[r["Type"]] = type_counts.get(r["Type"], 0) + 1
                breakdown = " | ".join(f"{k}: {v}" for k, v in sorted(type_counts.items()))
                print(f"           -> {breakdown}")
    print("="*80 + "\n")
if __name__ == "__main__":
    run_step1dc_conservancy_benchmark()
