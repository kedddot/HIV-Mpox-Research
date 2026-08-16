import os
import sys
import time
import csv
from datetime import datetime

# =============================================================================
# REPRESENTATIVE ALLELE FREQUENCY DATABASE
# Simulated 2-field notation frequencies (normally downloaded from IEDB)
# Missing alleles will safely default to 0.0 as per methodology.
# =============================================================================
MOCK_ALLELE_DB = {
    "HLA-A*02:01": 0.27, "HLA-A*24:02": 0.15, "HLA-A*01:01": 0.10,
    "HLA-B*07:02": 0.12, "HLA-B*08:01": 0.09, "HLA-B*35:01": 0.11,
    "HLA-C*07:01": 0.14, "HLA-C*04:01": 0.12,
    "HLA-DRB1*01:01": 0.09, "HLA-DRB1*15:01": 0.11, "HLA-DRB1*03:01": 0.08
}

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def calculate_population_coverage(alleles, freq_db):
    """
    Implements the precise mathematical formula for population coverage:
    1. Sum allele frequencies per locus (capped at 1.0)
    2. Per-locus coverage = 1 - (1 - s)^2
    3. Overall coverage = 1 - product of (1 - locus_coverage)
    """
    if not alleles:
        return 0.0
        
    # Group alleles by locus (e.g., "A", "B", "C", "DRB1")
    loci_sums = {}
    for allele in alleles:
        # Extract locus from "HLA-A*02:01" -> "A"
        try:
            locus = allele.split('-')[1].split('*')[0]
        except IndexError:
            continue
            
        freq = freq_db.get(allele, 0.0) # Missing treated as 0
        loci_sums[locus] = loci_sums.get(locus, 0.0) + freq

    # Calculate per-locus coverage
    locus_coverages = []
    for locus, s in loci_sums.items():
        s = min(1.0, s) # Cap at 1.0
        locus_cov = 1.0 - (1.0 - s)**2
        locus_coverages.append(locus_cov)

    # Combine loci assuming independence
    if not locus_coverages:
        return 0.0
        
    not_covered = 1.0
    for cov in locus_coverages:
        not_covered *= (1.0 - cov)
        
    overall_coverage = 1.0 - not_covered
    return overall_coverage * 100.0 # Return as percentage

def assign_mock_alleles(pep_type):
    """Assigns representative alleles for testing the mathematical model."""
    if pep_type == "MHC-I":
        return ["HLA-A*02:01", "HLA-A*24:02", "HLA-B*07:02", "HLA-C*07:01"]
    elif pep_type == "MHC-II":
        return ["HLA-DRB1*01:01", "HLA-DRB1*15:01", "HLA-DRB1*03:01"]
    return []

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def run_step1f_population_coverage():
    start_time = time.time()
    COVERAGE_THRESHOLD = 90.0
    
    # DYNAMIC PATHING
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1E", "Phase1Eb", "Filtered")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase1F")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)

    print("\n" + "="*80)
    print(f"{'PHASE 1F: POPULATION COVERAGE ANALYSIS (IEDB MATH MODEL)':^80}")
    print("="*80)

    try:
        csv_files = [f for f in os.listdir(input_folder) if f.endswith(".csv")]
        latest_csv = max([os.path.join(input_folder, f) for f in csv_files], key=os.path.getctime)
    except ValueError:
        print("[ERROR] No input CSV found in Phase 1Eb Filtered directory.")
        return

    print(f"[INFO] Source File        : {os.path.basename(latest_csv)}")
    print(f"[INFO] MHC-I Threshold    : Rank <= 2.0%")
    print(f"[INFO] MHC-II Threshold   : Rank <= 10.0%")
    print(f"[INFO] Coverage Threshold : >= 90.0%")
    print("-" * 80)

    raw_data = []
    filtered_data = []
    
    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        original_fields = reader.fieldnames
        # Ensure we don't duplicate field names if running multiple times
        base_fields = [fn for fn in original_fields if fn not in ["Projected_Coverage", "Coverage_Status", "Exclusion_Reason"]]
        fieldnames = base_fields + ["Projected_Coverage", "Coverage_Status", "Exclusion_Reason"]
        
        rows = list(reader)
        total_rows = len(rows)

        for i, row in enumerate(rows):
            pep_type = row.get('Type', 'Unknown')
            rank_str = row.get('Percentile_Rank') or row.get('Rank') or "999"
            try:
                rank = float(rank_str)
            except ValueError:
                rank = 999.0
            
            is_valid = True
            exclusion_reason = []

            # 1. Percentile Rank Validation (Skipped for B-cells as they use BepiPred)
            if pep_type == "MHC-I" and rank > 2.0:
                is_valid = False
                exclusion_reason.append(f"MHC-I Rank {rank}% > 2%")
            elif pep_type == "MHC-II" and rank > 10.0:
                is_valid = False
                exclusion_reason.append(f"MHC-II Rank {rank}% > 10%")

            # 2. Population Coverage Calculation (HLA-dependent)
            coverage_val = 100.0 # Default for non-HLA peptides (e.g. B-cell)
            if pep_type in ["MHC-I", "MHC-II"] and is_valid:
                alleles = assign_mock_alleles(pep_type)
                coverage_val = calculate_population_coverage(alleles, MOCK_ALLELE_DB)
                
                if coverage_val < COVERAGE_THRESHOLD:
                    is_valid = False
                    exclusion_reason.append(f"Coverage {coverage_val:.1f}% < 90%")
            
            # Progress Tracking
            elapsed = format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{total_rows:03d} | Calculating Coverage | Elapsed: {elapsed}")
            sys.stdout.flush()
            
            clean_row = {k: row.get(k, '') for k in base_fields}
            clean_row.update({
                "Projected_Coverage": f"{coverage_val:.2f}%" if pep_type in ["MHC-I", "MHC-II"] else "N/A (B-cell)",
                "Coverage_Status": "PASSED" if is_valid else "FAILED",
                "Exclusion_Reason": " | ".join(exclusion_reason) if not is_valid else "None"
            })
            
            raw_data.append(clean_row)
            if is_valid:
                filtered_data.append(clean_row)

    print(f"\n[INFO] Coverage analysis complete.")

    # EXPORT
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    if raw_data:
        with open(os.path.join(raw_dir, f"Phase1F_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(raw_data)
            
    if filtered_data:
        with open(os.path.join(filt_dir, f"Phase1F_Elite_Vaccine_Candidates_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(filtered_data)

    survivor_rate = (len(filtered_data) / len(raw_data)) * 100 if raw_data else 0

    print("-" * 80)
    print(f"[SUCCESS] Total Evaluated      : {len(raw_data)}")
    print(f"[SUCCESS] Elite Global Binders : {len(filtered_data)} ({survivor_rate:.1f}% Pass Rate)")
    print(f"[SUCCESS] Total Time           : {format_time(time.time() - start_time)}")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_step1f_population_coverage()