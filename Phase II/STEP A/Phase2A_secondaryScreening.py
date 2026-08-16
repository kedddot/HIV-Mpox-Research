import os
import csv
import sys
import time
from datetime import datetime
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =============================================================================
# ACADEMIC FORMATTING & BENCHMARKING UTILITIES (Retained from Original)
# =============================================================================

def print_banner(text):
    """Prints a formal academic section header."""
    print("\n" + "="*125)
    print(f"{text:^125}")
    print("="*125)

def format_time(seconds):
    """Formats execution time into a readable MM:SS string."""
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

# =============================================================================
# MOCK EXTERNAL PREDICTORS (Replace with actual API/Subprocess calls)
# =============================================================================
class ExternalPredictors:
    @staticmethod
    def run_toxinpred(sequence):
        return "non-toxin"

    @staticmethod
    def run_mhc_binding(junction_peptide, mhc_class):
        return {"rank": 5.0, "toxinpred": "non-toxin"}

    @staticmethod
    def run_allertop(sequence):
        return False

    @staticmethod
    def run_allergenfp(sequence):
        return False

    @staticmethod
    def run_solpro(sequence):
        return 0.85 

    @staticmethod
    def run_camsol(sequence):
        return 5

# =============================================================================
# BIOCHEMICAL & METHODOLOGICAL EVALUATION ENGINE
# =============================================================================

def get_sliding_windows(seq, window_size):
    return [seq[i:i+window_size] for i in range(len(seq) - window_size + 1)]

def check_hydrophobicity(seq):
    """Fails if any 8, 12, or 15 aa window is > 80% hydrophobic."""
    hydrophobic_residues = set(['A', 'C', 'F', 'I', 'L', 'M', 'V', 'W'])
    for w_size in [8, 12, 15]:
        for window in get_sliding_windows(seq, w_size):
            h_count = sum(1 for aa in window if aa in hydrophobic_residues)
            if (h_count / w_size) > 0.80:
                return False
    return True

def generate_junctions(seq, linkers=None):
    """
    Identifies predefined Phase I linkers in the sequence and extracts 
    a 20-amino-acid junction window centered on each linker.
    """
    if linkers is None:
        linkers = ["AAY", "GPGPG", "KK", "EAAAK"]
        
    junctions = []
    
    for linker in linkers:
        start_idx = 0
        while True:
            idx = seq.find(linker, start_idx)
            if idx == -1:
                break
                
            linker_len = len(linker)
            flank_needed = 20 - linker_len
            left_flank = flank_needed // 2
            right_flank = flank_needed - left_flank
            
            window_start = max(0, idx - left_flank)
            window_end = min(len(seq), idx + linker_len + right_flank)
            
            junction_peptide = seq[window_start:window_end]
            
            if len(junction_peptide) >= 8:
                junctions.append(junction_peptide)
                
            start_idx = idx + 1
            
    return list(set(junctions))

def evaluate_candidate(filename, seq):
    """Executes the strict Phase II methodology checks on a given sequence."""
    ana = ProteinAnalysis(seq)
    instability_idx = ana.instability_index()
    gravy_val = ana.gravy()
    pI = ana.isoelectric_point()
    
    result = {
        "Variant": filename,
        "Len": len(seq),
        "STAB_IDX": round(instability_idx, 2),
        "GRAVY": round(gravy_val, 4),
        "pI": round(pI, 2),
        "Viable": "NO",
        "Rejection_Reasons": [],
        "Review_Flags": []
    }
    
    # 1. TOXICITY SCREENING
    if not check_hydrophobicity(seq):
        result["Rejection_Reasons"].append("Hydrophobicity >80% in window")
        
    junctions = generate_junctions(seq)
    for j in junctions:
        if 8 <= len(j) <= 20 and j.count('C') > 2:
            result["Rejection_Reasons"].append("Junction Cys > 2")
            break
            
    if ExternalPredictors.run_toxinpred(seq) == "toxin":
        result["Rejection_Reasons"].append("ToxinPred: Toxin")
        
    for j in junctions:
        mhc_i = ExternalPredictors.run_mhc_binding(j, "MHC_I")
        mhc_ii = ExternalPredictors.run_mhc_binding(j, "MHC_II")
        if (mhc_i["rank"] <= 2.0 or mhc_ii["rank"] <= 10.0) and mhc_i["toxinpred"] == "toxin":
            result["Rejection_Reasons"].append("Toxic Junction (MHC)")
            break

    # 2. ALLERGENICITY SCREENING
    # NOTE: Per methodology, Surface_Charge only "flags the sequence for review" and is
    # NOT part of the final consensus rejection rule. The final rejection rule is based
    # solely on: (a) both AllerTop + AllergenFP positive, or (b) QN_Ratio > 0.30 combined
    # with >=1 predictor positive, or (c) a junction positive on both predictors.
    qn_ratio = (seq.count('Q') + seq.count('N')) / len(seq)
    surface_charge = sum(seq.count(aa) for aa in ['D', 'E', 'H', 'K'])
    charge_threshold = max(4, round(len(seq) * 0.02))
    
    if surface_charge > charge_threshold:
        result["Review_Flags"].append(
            f"High Surface Charge ({surface_charge} residues > threshold {charge_threshold})"
        )
        
    aller_top = ExternalPredictors.run_allertop(seq)
    aller_fp = ExternalPredictors.run_allergenfp(seq)
    allergen_count = sum([aller_top, aller_fp])
    
    if allergen_count == 2:
        result["Rejection_Reasons"].append("Allergenic (Both Predictors)")
    elif qn_ratio > 0.30 and allergen_count >= 1:
        result["Rejection_Reasons"].append("Allergenic (QN > 0.30 + 1 Predictor)")
        
    for j in junctions:
        if ExternalPredictors.run_allertop(j) and ExternalPredictors.run_allergenfp(j):
            result["Rejection_Reasons"].append("Allergenic Junction")
            break

    # Junctional peptide surface charge: raw Surface_Charge > 4 is a review flag only
    # ("automatically deprioritized", not an automatic reject) per methodology.
    for j in junctions:
        j_charge = sum(j.count(aa) for aa in ['D', 'E', 'H', 'K'])
        if j_charge > 4:
            result["Review_Flags"].append(f"Junction High Surface Charge ({j_charge} in '{j}')")

    # 3. SOLUBILITY & STABILITY SCREENING
    if instability_idx > 40:
        result["Rejection_Reasons"].append("Instability > 40")
    if gravy_val > 0.4:
        result["Rejection_Reasons"].append("GRAVY > 0.4")
    if 6.5 <= pI <= 7.5:
        result["Rejection_Reasons"].append("pI near 7.0")
        
    solpro_score = ExternalPredictors.run_solpro(seq)
    camsol_patch = ExternalPredictors.run_camsol(seq)
    
    if solpro_score < 0.5:
        result["Rejection_Reasons"].append("SOLpro Insoluble")
    if camsol_patch > 8 and solpro_score < 0.5:
        result["Rejection_Reasons"].append("CamSol Aggregation > 8aa + Insoluble")

    # FINAL VERDICT (Review_Flags never factor into Viable status)
    if not result["Rejection_Reasons"]:
        result["Viable"] = "YES"
        result["Rejection_Reasons"] = "None (Optimal)"
    else:
        result["Rejection_Reasons"] = " | ".join(result["Rejection_Reasons"])

    result["Review_Flags"] = " | ".join(result["Review_Flags"]) if result["Review_Flags"] else "None"

    return result

# =============================================================================
# MAIN EXECUTION THREAD
# =============================================================================

def run_step2a_comprehensive_screening():
    start_time = time.time()

    # 1. DYNAMIC DIRECTORY RESOLUTION
    # NOTE: A fixed "../../.." traversal is fragile -- it silently breaks
    # whenever a script moves to a folder at a different nesting depth.
    # Instead, walk upward from the script's own location until we find
    # the "Research" anchor folder, matching the logic used in Step 2B.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = script_dir
    while os.path.basename(project_root) != "Research":
        parent = os.path.dirname(project_root)
        if parent == project_root:
            print(f"\n[FATAL ERROR] Could not locate a 'Research' anchor folder above: {script_dir}")
            sys.exit(1)
        project_root = parent

    input_fasta_path = "/Users/nek/Desktop/School/Research/Phase 1/Step_Outputs/Phase1G/Phase1G_Constructs_2026-07-22_2251.fasta"
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA")
    
    raw_out_dir = os.path.join(output_base, "Raw")
    filt_out_dir = os.path.join(output_base, "Filtered")
    
    os.makedirs(raw_out_dir, exist_ok=True)
    os.makedirs(filt_out_dir, exist_ok=True)

    if not os.path.exists(input_fasta_path):
        print(f"\n[FATAL ERROR] Variant FASTA file not found: {input_fasta_path}")
        sys.exit(1)

    # 2. EXPERIMENTAL INITIALIZATION (Retained & Updated)
    print_banner("PHASE 2 STEP A: COMPREHENSIVE SECONDARY SCREENING DASHBOARD")
    print(f"[INFO] Initialization Time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Input Source        : {input_fasta_path}")
    print(f"[INFO] Screening Standards : Rigorous Phase II Methodology Integrated")
    print("-" * 125)

    raw_results = []

    # 3. ANALYTICAL EXECUTION
    # Parse all records out of the single multi-sequence FASTA file.
    variants = []
    current_name = None
    current_seq_lines = []

    with open(input_fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    variants.append((current_name, "".join(current_seq_lines).upper()))
                current_name = line[1:].strip()
                current_seq_lines = []
            else:
                current_seq_lines.append(line)
        if current_name is not None:
            variants.append((current_name, "".join(current_seq_lines).upper()))

    if not variants:
        print("[WARNING] No FASTA variants found in Phase 1G. Please check Step 1G execution.")
        return

    for filename, seq in variants:
        candidate_data = evaluate_candidate(filename, seq)
        raw_results.append(candidate_data)

    # 4. TERMINAL OUTPUT: RAW DATA TABLE
    print(f"{'VARIANT':<22} | {'LEN':<4} | {'STAB':<6} | {'GRAVY':<6} | {'pI':<5} | {'VIABLE':<6} | {'REASONS':<48} | {'REVIEW_FLAGS'}")
    print("-" * 125)
    for r in raw_results:
        # Truncate reasons for terminal display to avoid massive text wrap
        short_reasons = (r['Rejection_Reasons'][:45] + '...') if len(r['Rejection_Reasons']) > 45 else r['Rejection_Reasons']
        short_flags = (r['Review_Flags'][:35] + '...') if len(r['Review_Flags']) > 35 else r['Review_Flags']
        print(f"{r['Variant'][:22]:<22} | {r['Len']:<4} | {r['STAB_IDX']:<6} | {r['GRAVY']:<6} | {r['pI']:<5} | {r['Viable']:<6} | {short_reasons:<48} | {short_flags}")

    # 5. DATA FILTERING & RANKING
    filtered_list = [r for r in raw_results if r['Viable'] == "YES"]
    filtered_list.sort(key=lambda x: x['STAB_IDX'])

    # 6. TERMINAL OUTPUT: FILTERED & RANKED TABLE
    print_banner("STEP 2A: FILTERED OPTIMAL CANDIDATES (RANKED BY STABILITY INDEX)")
    if filtered_list:
        print(f"{'RANK':<5} | {'VARIANT':<22} | {'STABILITY INDEX':<18} | {'GRAVY':<10} | {'STATUS'}")
        print("-" * 125)
        for i, f in enumerate(filtered_list):
            status_label = "[OPTIMAL]" if i == 0 else "[ACCEPTED]"
            print(f"{i+1:<5} | {f['Variant']:<22} | {f['STAB_IDX']:<18} | {f['GRAVY']:<10} | {status_label}")
    else:
        print(f"{'--- NO VARIANTS PASSED THE COMPREHENSIVE SCREENING CRITERIA ---':^125}")

    # 7. DATA EXPORT (Dual CSV Logging) (Retained)
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    raw_path = os.path.join(raw_out_dir, f"Step2A_Raw_Physicochemical_{ts}.csv")
    filt_path = os.path.join(filt_out_dir, f"Step2A_Filtered_Ranked_{ts}.csv")

    with open(raw_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=raw_results[0].keys())
        writer.writeheader()
        writer.writerows(raw_results)

    if filtered_list:
        with open(filt_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=filtered_list[0].keys())
            writer.writeheader()
            writer.writerows(filtered_list)

    # 8. FINAL BENCHMARK SUMMARY (Retained)
    exec_time = format_time(time.time() - start_time)
    print_banner("PHASE 2 STEP A COMPLETE")
    print(f"[SUCCESS] Total Analyzed      : {len(raw_results)} variants")
    print(f"[SUCCESS] Viable Candidates   : {len(filtered_list)}")
    print(f"[SUCCESS] Execution Time      : {exec_time}")
    print(f"[INFO] Raw Log: {os.path.relpath(raw_path, project_root)}")
    print(f"[INFO] Filtered Log: {os.path.relpath(filt_path, project_root)}")
    print("="*125 + "\n")

if __name__ == "__main__":
    run_step2a_comprehensive_screening()
