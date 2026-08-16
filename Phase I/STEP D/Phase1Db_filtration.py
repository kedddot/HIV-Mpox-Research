import os, sys, time, re, csv, requests
from datetime import datetime
from collections import defaultdict

# =============================================================================
# METHODOLOGY / ALLELE SCOPE NOTE
# =============================================================================
MHCI_ALLELE = "HLA-A*02:01"
MHCII_ALLELE = "HLA-DRB1*01:01"

MHCI_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
MHCII_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"
BCELL_URL = "https://tools-cluster-interface.iedb.org/tools_api/bcell/"

def print_banner(text): print(f"\n{'='*80}\n{text:^80}\n{'='*80}")

def get_gravy(pep):
    hydro = {'A': 1.8, 'L': 3.8, 'I': 4.5, 'V': 4.2, 'F': 2.8, 'M': 1.9, 'C': 2.5, 'G': -0.4,
             'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'N': -3.5, 'Q': -3.5,
             'D': -3.5, 'E': -3.5, 'K': -3.9, 'R': -4.5, 'H': -3.2}
    return sum(hydro.get(aa, 0) for aa in pep) / len(pep)

def classify_bcell_tier(mean_score, pct_above):
    if mean_score >= 0.60 and pct_above >= 75.0: return "High"
    elif mean_score >= 0.50 and pct_above >= 50.0: return "Medium"
    elif mean_score >= 0.45 and pct_above >= 37.5: return "Deprioritized"
    return "Excluded"

def find_column(header, keywords):
    for kw in keywords:
        for i, col in enumerate(header):
            if kw.lower() in col.lower(): return i
    return None

# --- NEW: Live Micro-Update Tracker ---
def print_status(current, total, target, start_time, kept, action):
    elapsed = time.time() - start_time
    # Pad with spaces to overwrite any leftover characters from previous longer lines
    msg = f"[PROCESS] {current:02d}/{total:02d} | Target: {target:<12} | Action: {action:<20} | Elapsed: {elapsed:5.1f}s | Kept: {kept}"
    sys.stdout.write("\r" + " " * 100)
    sys.stdout.write(f"\r{msg}")
    sys.stdout.flush()

def run_step1db_optimized():
    start_time = time.time()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    fasta_folder = os.path.join(project_root, "Step_Outputs", "Phase1C", "Filtered_Antigenicity")
    identification_folder = os.path.join(project_root, "Step_Outputs", "Phase1D", "Phase1Da")
    output_folder = os.path.join(project_root, "Step_Outputs", "Phase1D", "Phase1Db")
    os.makedirs(output_folder, exist_ok=True)

    print_banner("PHASE 1Db: IEDB AFFINITY & SOLUBILITY FILTER")

    if not os.path.exists(identification_folder):
        print(f"\n[ERROR] Phase 1Da directory not found at: {identification_folder}")
        return
    da_files = [f for f in os.listdir(identification_folder) if f.endswith(".csv")]
    if not da_files:
        print(f"\n[ERROR] No Phase 1Da candidate library (.csv) found in: {identification_folder}")
        return
    latest_da = os.path.join(identification_folder, sorted(da_files)[-1])
    
    print(f"[INFO] Candidate Library    : {os.path.basename(latest_da)}")
    print(f"[INFO] MHC-I Allele         : {MHCI_ALLELE}")
    print(f"[INFO] MHC-II Allele        : {MHCII_ALLELE}")
    print("-" * 80)

    candidates_by_variant = defaultdict(list)
    with open(latest_da, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader: candidates_by_variant[row["Variant"]].append(row)

    fasta_files = sorted([f for f in os.listdir(fasta_folder) if f.endswith(".fasta")])

    all_results = []
    skipped = {"MHC-I": 0, "MHC-II": 0, "B-cell": 0}
    api_errors = 0
    
    # NEW: Using a Requests Session to speed up server handshakes
    session = requests.Session()

    for i, f_name in enumerate(fasta_files):
        target = f_name.split('_Var')[0]
        my_candidates = candidates_by_variant.get(f_name, [])
        if not my_candidates: continue

        with open(os.path.join(fasta_folder, f_name), "r") as f:
            seq = "".join([l.strip() for l in f if not l.startswith(">")])
            clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())
        if len(clean_seq) < 9: continue

        bcell_candidates = [c for c in my_candidates if c["Type"] == "B-cell"]

        # ---------------- MHC-I: 9-mer and 10-mer ----------------
        for length in (9, 10):
            print_status(i+1, len(fasta_files), target, start_time, len(all_results), f"Query MHC-I ({length}-mer)...")
            payload = {'method': 'recommended', 'sequence_text': clean_seq, 'allele': MHCI_ALLELE, 'length': str(length)}
            try:
                response = session.post(MHCI_URL, data=payload, timeout=90)
                if response.status_code == 200:
                    lines = response.text.strip().split('\n')
                    if len(lines) > 1:
                        header = lines[0].split('\t')
                        r_idx = find_column(header, ['percentile_rank', 'rank'])
                        p_idx = find_column(header, ['peptide'])
                        if r_idx is not None and p_idx is not None:
                            for line in lines[1:]:
                                cols = line.split('\t')
                                pep = cols[p_idx]
                                try: rank = float(cols[r_idx])
                                except: continue
                                
                                if rank <= 1.0:
                                    gravy = get_gravy(pep)
                                    if gravy < 0.2:
                                        all_results.append({"Target": target, "Variant": f_name, "Type": "MHC-I", "Length": length, "Peptide": pep, "GRAVY": round(gravy, 3), "Percentile_Rank": rank, "mean_BepiPred": "", "pct_above": "", "Bcell_Tier": ""})
                                    else: skipped["MHC-I"] += 1
                else: api_errors += 1
                time.sleep(1.0) 
            except Exception as e:
                api_errors += 1
                print(f"\n[WARN] MHC-I request failed for {f_name} (len {length}): Timeout/Server Error.")
                continue

        # ---------------- MHC-II: 15-mer ----------------
        print_status(i+1, len(fasta_files), target, start_time, len(all_results), "Query MHC-II...")
        payload = {'method': 'recommended', 'sequence_text': clean_seq, 'allele': MHCII_ALLELE, 'length': '15'}
        try:
            response = session.post(MHCII_URL, data=payload, timeout=90)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                if len(lines) > 1:
                    header = lines[0].split('\t')
                    r_idx = find_column(header, ['percentile_rank', 'adjusted_rank', 'rank'])
                    p_idx = find_column(header, ['peptide'])
                    if r_idx is not None and p_idx is not None:
                        for line in lines[1:]:
                            cols = line.split('\t')
                            pep = cols[p_idx]
                            try: rank = float(cols[r_idx])
                            except: continue
                            
                            if rank <= 10.0:
                                gravy = get_gravy(pep)
                                if gravy < 0.2:
                                    all_results.append({"Target": target, "Variant": f_name, "Type": "MHC-II", "Length": 15, "Peptide": pep, "GRAVY": round(gravy, 3), "Percentile_Rank": rank, "mean_BepiPred": "", "pct_above": "", "Bcell_Tier": ""})
                                else: skipped["MHC-II"] += 1
            else: api_errors += 1
            time.sleep(1.0)
        except Exception as e:
            api_errors += 1
            print(f"\n[WARN] MHC-II request failed for {f_name}: Timeout/Server Error.")

        # ---------------- B-cell: 16-mer ----------------
        if bcell_candidates:
            print_status(i+1, len(fasta_files), target, start_time, len(all_results), "Query B-cell...")
            payload = {'method': 'Bepipred-2.0', 'sequence_text': clean_seq}
            try:
                response = session.post(BCELL_URL, data=payload, timeout=90)
                if response.status_code == 200:
                    lines = response.text.strip().split('\n')
                    if len(lines) > 1:
                        header = lines[0].split('\t')
                        score_idx = find_column(header, ['score'])
                        if score_idx is not None:
                            residue_scores = []
                            for line in lines[1:]:
                                cols = line.split('\t')
                                try: residue_scores.append(float(cols[score_idx]))
                                except: continue
                            for cand in bcell_candidates:
                                start = int(cand["Start_Position"])
                                window = residue_scores[start:start + 16]
                                if len(window) < 16: continue
                                mean_score = sum(window) / len(window)
                                pct_above = 100.0 * sum(1 for s in window if s >= 0.50) / len(window)
                                tier = classify_bcell_tier(mean_score, pct_above)
                                
                                if tier in ("High", "Medium"):
                                    gravy = get_gravy(cand["Peptide"])
                                    if gravy < 0.2:
                                        all_results.append({"Target": target, "Variant": f_name, "Type": "B-cell", "Length": 16, "Peptide": cand["Peptide"], "GRAVY": round(gravy, 3), "Percentile_Rank": "", "mean_BepiPred": round(mean_score, 3), "pct_above": round(pct_above, 2), "Bcell_Tier": tier})
                                    else: skipped["B-cell"] += 1
                                else: skipped["B-cell"] += 1
                else: api_errors += 1
                time.sleep(1.0)
            except Exception as e:
                api_errors += 1
                print(f"\n[WARN] B-cell request failed for {f_name}: Timeout/Server Error.")

    print()  # newline after progress bar

    # 5. Export
    out_file = os.path.join(output_folder, f"Phase1Db_Elite_{datetime.now().strftime('%Y%m%d')}.csv")
    fieldnames = ["Target", "Variant", "Type", "Length", "Peptide", "GRAVY",
                  "Percentile_Rank", "mean_BepiPred", "pct_above", "Bcell_Tier"]
    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    print_banner("FILTRATION COMPLETE")
    n_mhci = sum(1 for r in all_results if r["Type"] == "MHC-I")
    n_mhcii = sum(1 for r in all_results if r["Type"] == "MHC-II")
    n_bcell = sum(1 for r in all_results if r["Type"] == "B-cell")
    print(f"[SUCCESS] {len(all_results)} total high-affinity/high-tier soluble epitopes saved.")
    print(f"          - MHC-I  : {n_mhci:<5} (dropped -- rank/GRAVY: {skipped['MHC-I']})")
    print(f"          - MHC-II : {n_mhcii:<5} (dropped -- rank/GRAVY: {skipped['MHC-II']})")
    print(f"          - B-cell : {n_bcell:<5} (dropped -- tier/GRAVY: {skipped['B-cell']})")
    if api_errors:
        print(f"[WARN] {api_errors} API request(s) failed or returned a non-200 status -- see warnings above.")
    print(f"[INFO] Log file: {os.path.basename(out_file)}\n")

if __name__ == "__main__":
    run_step1db_optimized()