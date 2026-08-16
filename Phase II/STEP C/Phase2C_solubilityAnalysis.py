import os
import csv
import sys
import time
from datetime import datetime

# =============================================================================
# ACADEMIC FORMATTING UTILITIES
# =============================================================================

def print_banner(text):
    print("\n" + "="*110)
    print(f"{text:^110}")
    print("="*110)

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

# =============================================================================
# MOCK EXTERNAL 3D STRUCTURAL TOOLS (Replace with PDB parsers / server APIs)
# =============================================================================

class MockStructuralTools:
    @staticmethod
    def run_solpro(sequence):
        # Mock: Returns SOLpro overall predicted solubility probability
        return 0.48 

    @staticmethod
    def analyze_aggrescan3d(pdb_file):
        # Mock: Returns max contiguous patch length and if it is in a high pLDDT region
        return {"max_patch_length": 10, "in_high_pLDDT_region": True}

    @staticmethod
    def analyze_sasa(pdb_file):
        # Mock: Returns hydrophobic surface fraction and max exposed patch area (Å²)
        return {"hydrophobic_fraction": 0.28, "exposed_patch_area": 260.0}

    @staticmethod
    def run_foldx(pdb_file):
        # Mock: Returns destabilizing ΔΔG
        return 1.8 

    @staticmethod
    def generate_apbs_map(pdb_file, output_dir):
        # Mock: Generates electrostatic map
        map_path = os.path.join(output_dir, "electrostatic_surface.dx")
        with open(map_path, 'w') as f: f.write("Mock APBS Data")
        return map_path

# =============================================================================
# SOLUBILITY & STRUCTURAL INTEGRITY ENGINE
# =============================================================================

def evaluate_structural_solubility(sequence, pdb_file, output_dir):
    """Evaluates the candidate strictly against Phase II Step C decision rules."""
    
    # 1. Gather sequence and structural metrics
    solpro_prob = MockStructuralTools.run_solpro(sequence)
    agg_data = MockStructuralTools.analyze_aggrescan3d(pdb_file)
    sasa_data = MockStructuralTools.analyze_sasa(pdb_file)
    ddg = MockStructuralTools.run_foldx(pdb_file)
    
    # Mock generating supplementary files
    MockStructuralTools.generate_apbs_map(pdb_file, output_dir)
    
    result = {
        "SOLpro_Prob": solpro_prob,
        "Agg_Patch_Len": agg_data["max_patch_length"],
        "High_pLDDT_Patch": agg_data["in_high_pLDDT_region"],
        "Hydrophobic_Fraction": sasa_data["hydrophobic_fraction"],
        "Exposed_Area": sasa_data["exposed_patch_area"],
        "FoldX_ddG": ddg,
        "Status": "PASS",
        "Reason_Code": []
    }

    # 2. Decision Logic: Hard Rejection
    is_insoluble = solpro_prob < 0.50
    has_confident_patch = agg_data["max_patch_length"] > 8 and agg_data["in_high_pLDDT_region"]
    
    if is_insoluble and has_confident_patch:
        result["Status"] = "REJECT"
        result["Reason_Code"].append("SOLpro <0.5 & Confident Aggregation Patch >8aa")
        
    # 3. Decision Logic: Flag for Review (if not already rejected)
    if result["Status"] != "REJECT":
        if sasa_data["hydrophobic_fraction"] > 0.25:
            result["Status"] = "REVIEW"
            result["Reason_Code"].append("Hydrophobic Fraction > 0.25")
        if sasa_data["exposed_patch_area"] > 250.0:
            result["Status"] = "REVIEW"
            result["Reason_Code"].append("Exposed Hydrophobic Area > 250 Å²")
        if ddg > 1.5:
            result["Status"] = "REVIEW"
            result["Reason_Code"].append("FoldX ΔΔG > +1.5 kcal/mol")
            
    if not result["Reason_Code"]:
        result["Reason_Code"].append("All structural parameters optimal")

    result["Reason_Code"] = " | ".join(result["Reason_Code"])
    return result

def run_step2c_solubility_analysis():
    start_time = time.time()

    # 1. DYNAMIC PATHING & ORGANIZATION
    # NOTE: A fixed "../../.." traversal is fragile -- it silently breaks
    # whenever a script moves to a folder at a different nesting depth.
    # Instead, walk upward from the script's own location until we find
    # the "Research" anchor folder, matching Steps 2A and 2B.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = script_dir
    while os.path.basename(project_root) != "Research":
        parent = os.path.dirname(project_root)
        if parent == project_root:
            print(f"\n[FATAL ERROR] Could not locate a 'Research' anchor folder above: {script_dir}")
            sys.exit(1)
        project_root = parent

    input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    # NOTE: Phase 1G no longer outputs one FASTA file per variant into a
    # "Variants" folder -- it outputs a single multi-sequence FASTA file,
    # matching what Steps 2A and 2B read from.
    variant_fasta_path = os.path.join(project_root, "Phase 1", "Step_Outputs", "Phase1G", "Phase1G_Constructs_2026-07-22_2251.fasta")
    pdb_input_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB", "Tertiary_Structure")
    
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC")
    archive_dir = os.path.join(output_base, "Supplementary_Archive")
    
    os.makedirs(output_base, exist_ok=True)
    os.makedirs(archive_dir, exist_ok=True)

    # 2. IDENTIFY TARGET
    print(f"[INFO] Resolved Project Root : {project_root}")
    print(f"[INFO] Looking for Filtered CSVs in : {input_csv_dir}")

    if not os.path.isdir(input_csv_dir):
        print(f"[ERROR] Filtered output folder does not exist: {input_csv_dir}")
        print("[ERROR] Run Step 2A first so this folder gets created.")
        return

    csv_files = sorted([f for f in os.listdir(input_csv_dir) if f.endswith(".csv")])
    if not csv_files:
        print(f"[ERROR] No filtered candidate CSVs found in: {input_csv_dir}")
        print("[ERROR] Step 2A ran but produced zero viable candidates -- check Rejection_Reasons in the Raw log.")
        return
    
    with open(os.path.join(input_csv_dir, csv_files[-1]), 'r') as f:
        reader = list(csv.DictReader(f))
        winner_name = reader[0]['Variant'] 

    print_banner("PHASE 2 STEP C: 3D SOLUBILITY & STRUCTURAL INTEGRITY")
    print(f"[INFO] Target Variant : {winner_name}")
    print(f"[INFO] Methodology    : PDB2PQR, APBS, Aggrescan3D, CamSol, FoldX")
    print("-" * 110)

    # 3. LOAD SEQUENCE & MOCK PDB
    # Parse all records out of the single multi-sequence FASTA file and
    # match by header (exact "Variant" text from Step 2A's CSV first,
    # then falling back to the name before any "|" metadata annotation).
    if not os.path.isfile(variant_fasta_path):
        print(f"[ERROR] Variant FASTA file not found: {variant_fasta_path}")
        return

    variants = {}
    current_name = None
    current_seq_lines = []
    with open(variant_fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    variants[current_name] = "".join(current_seq_lines).upper()
                current_name = line[1:].strip()
                current_seq_lines = []
            else:
                current_seq_lines.append(line)
        if current_name is not None:
            variants[current_name] = "".join(current_seq_lines).upper()

    clean_name = winner_name.split("|")[0].strip()
    sequence = variants.get(winner_name) or variants.get(clean_name)

    if sequence is None:
        print(f"[ERROR] Could not find a matching header for '{winner_name}' in: {variant_fasta_path}")
        available = list(variants.keys())
        print(f"[ERROR] Headers actually present ({len(available)}): {available[:10]}{' ...' if len(available) > 10 else ''}")
        return

    safe_name = winner_name.split("|")[0].strip().replace(" ", "_")
    mock_pdb_path = os.path.join(pdb_input_dir, f"AF3_Target_{safe_name}.pdb")
    
    # 4. EXECUTE ANALYSIS
    results = evaluate_structural_solubility(sequence, mock_pdb_path, archive_dir)

    # 5. TERMINAL DASHBOARD
    print(f"{'METRIC':<25} | {'VALUE':<15} | {'THRESHOLD / TARGET'}")
    print("-" * 110)
    print(f"{'SOLpro Probability':<25} | {results['SOLpro_Prob']:<15.2f} | > 0.50")
    print(f"{'Max Aggregation Patch':<25} | {results['Agg_Patch_Len']:<15} | < 8 aa (pLDDT >= 70)")
    print(f"{'Hydrophobic Fraction':<25} | {results['Hydrophobic_Fraction']:<15.2f} | < 0.25")
    print(f"{'Exposed Hydro Area':<25} | {results['Exposed_Area']:<15.1f} | < 250 Å²")
    print(f"{'FoldX ΔΔG':<25} | {results['FoldX_ddG']:<15.2f} | < +1.5 kcal/mol")
    print("-" * 110)
    print(f"FINAL DECISION : [{results['Status']}] - {results['Reason_Code']}")
    print("-" * 110)

    # 6. DATA EXPORT & ARCHIVING
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    report_path = os.path.join(output_base, f"Step2C_Solubility_Report_{ts}.csv")
    
    csv_data = {"Variant": winner_name, **results}
    
    with open(report_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=csv_data.keys())
        w.writeheader()
        w.writerow(csv_data)

    total_time = format_time(time.time() - start_time)
    print_banner("STEP 2C COMPLETE")
    print(f"[SUCCESS] Structural Solubilty & Aggregation Propensity Checked.")
    print(f"[SUCCESS] Execution Time : {total_time}")
    print(f"[INFO] Report Saved      : {os.path.relpath(report_path, project_root)}")
    print(f"[INFO] PDBs & Maps Saved : {os.path.relpath(archive_dir, project_root)}")
    print("="*110 + "\n")

if __name__ == "__main__":
    run_step2c_solubility_analysis()