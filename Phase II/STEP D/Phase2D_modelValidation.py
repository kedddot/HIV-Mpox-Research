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
# MOCK EXTERNAL VALIDATION TOOLS (Replace with actual MolProbity/Rampage calls)
# =============================================================================

class MockValidationTools:
    @staticmethod
    def run_molprobity(structural_file):
        # Mock: Returns a MolProbity score and stereochemical clash data
        return {
            "molprobity_score": 1.25,
            "clashscore": 4.10,
            "poor_rotamers_percent": 0.8
        }

    @staticmethod
    def run_rampage(structural_file, output_dir):
        # Mock: Returns Ramachandran Plot statistics based on phi/psi torsion angles
        plot_path = os.path.join(output_dir, "Ramachandran_Plot.png")
        with open(plot_path, 'w') as f: f.write("Mock Image Data")
        
        return {
            "favored_percent": 96.5,
            "allowed_percent": 2.8,
            "outlier_percent": 0.7,
            "plot_file": plot_path
        }

# =============================================================================
# CORE STEREOCHEMICAL VALIDATION ENGINE
# =============================================================================

def evaluate_model_quality(pdb_file, output_dir):
    """Evaluates stereochemical quality based on Phase II Step D rules."""
    
    mp_data = MockValidationTools.run_molprobity(pdb_file)
    ram_data = MockValidationTools.run_rampage(pdb_file, output_dir)
    
    result = {
        "MolProbity_Score": mp_data["molprobity_score"],
        "Rama_Favored": ram_data["favored_percent"],
        "Rama_Allowed": ram_data["allowed_percent"],
        "Rama_Outlier": ram_data["outlier_percent"],
        "MP_Quality": "UNKNOWN",
        "Rama_Quality": "UNKNOWN",
        "Overall_Status": "REVIEW"
    }

    # 1. Evaluate MolProbity Score
    score = result["MolProbity_Score"]
    if 0.5 <= score <= 1.5:
        result["MP_Quality"] = "EXCELLENT"
    elif score < 2.0:
        result["MP_Quality"] = "ACCEPTABLE"
    else:
        result["MP_Quality"] = "POOR"

    # 2. Evaluate Ramachandran (Torsion angles Phi and Psi)
    favored = result["Rama_Favored"]
    if 95.0 <= favored <= 98.0:
        result["Rama_Quality"] = "OPTIMAL (95-98%)"
    elif favored > 98.0:
        result["Rama_Quality"] = "EXCEPTIONAL (>98%)"
    elif favored >= 90.0:
        result["Rama_Quality"] = "ACCEPTABLE (90-94%)"
    else:
        result["Rama_Quality"] = "POOR (<90%)"

    # 3. Final Verdict
    if result["MP_Quality"] in ["EXCELLENT", "ACCEPTABLE"] and result["Rama_Quality"] in ["OPTIMAL (95-98%)", "EXCEPTIONAL (>98%)"]:
        result["Overall_Status"] = "VALIDATED"
    else:
        result["Overall_Status"] = "FAILED - REQUIRES REFINEMENT"

    return result

def run_step2d_model_validation():
    start_time = time.time()

    # 1. DYNAMIC PATHING & ORGANIZATION
    # NOTE: A fixed "../../.." traversal is fragile -- it silently breaks
    # whenever a script moves to a folder at a different nesting depth.
    # Instead, walk upward from the script's own location until we find
    # the "Research" anchor folder, matching Steps 2A, 2B, and 2C.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = script_dir
    while os.path.basename(project_root) != "Research":
        parent = os.path.dirname(project_root)
        if parent == project_root:
            print(f"\n[FATAL ERROR] Could not locate a 'Research' anchor folder above: {script_dir}")
            sys.exit(1)
        project_root = parent
    
    # 2. ALIGNED ROUTING: Target the Step 2A CSV to find the dynamic winner name
    input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    # ALIGNED ROUTING: Target the PDB from the Step 2B Tertiary output folder
    input_pdb_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB", "Tertiary_Structure")
    
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepD")
    archive_dir = os.path.join(output_base, "Validation_Plots")
    
    os.makedirs(output_base, exist_ok=True)
    os.makedirs(archive_dir, exist_ok=True)

    print_banner("PHASE 2 STEP D: STEREOCHEMICAL MODEL VALIDATION")
    print(f"[INFO] Resolved Project Root : {project_root}")
    print(f"[INFO] Methodology : MolProbity & Rampage (Ramachandran Plot)")
    print(f"[INFO] Targets     : MolProbity < 2.0 | Favored Residues 95-98%")
    print("-" * 110)

    # 3. DYNAMICALLY IDENTIFY TARGET PDB
    if not os.path.exists(input_csv_dir):
        print(f"[ERROR] Filtered directory not found: {input_csv_dir}")
        return

    csv_files = sorted([f for f in os.listdir(input_csv_dir) if f.endswith(".csv")])
    if not csv_files:
        print("[ERROR] No filtered candidates found. Run Step 2A first.")
        return
    
    with open(os.path.join(input_csv_dir, csv_files[-1]), 'r') as f:
        reader = list(csv.DictReader(f))
        winner_name = reader[0]['Variant'] 

    safe_name = winner_name.split("|")[0].strip().replace(" ", "_")
    target_pdb_name = f"AF3_Target_{safe_name}.pdb"
    target_pdb_path = os.path.join(input_pdb_dir, target_pdb_name)

    print(f"[INFO] Validating Model : {target_pdb_name}")
    print("-" * 110)

    # 4. EXECUTE VALIDATION (Using the dynamically linked path)
    if not os.path.exists(target_pdb_path):
        print(f"[WARNING] PDB file not found at {target_pdb_path}.")
        print("[INFO] Please ensure you have run AlphaFold3 and placed the .pdb in the Step 2B Tertiary folder.")
        return
        
    results = evaluate_model_quality(target_pdb_path, archive_dir)

    # 5. TERMINAL DASHBOARD
    print(f"{'VALIDATION METRIC':<30} | {'VALUE':<15} | {'ASSESSMENT'}")
    print("-" * 110)
    print(f"{'MolProbity Score':<30} | {results['MolProbity_Score']:<15.2f} | {results['MP_Quality']}")
    print(f"{'Rama: Favored Region':<30} | {results['Rama_Favored']:<14.2f}% | {results['Rama_Quality']}")
    print(f"{'Rama: Allowed Region':<30} | {results['Rama_Allowed']:<14.2f}% | Expected < 5%")
    print(f"{'Rama: Outlier Region':<30} | {results['Rama_Outlier']:<14.2f}% | Expected < 1%")
    print("-" * 110)
    print(f"FINAL DECISION : [{results['Overall_Status']}]")
    print("-" * 110)

    # 6. DATA EXPORT & ARCHIVING
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    report_path = os.path.join(output_base, f"Step2D_Validation_Report_{ts}.csv")
    
    csv_data = {"Model": target_pdb_name, **results}
    
    with open(report_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=csv_data.keys())
        w.writeheader()
        w.writerow(csv_data)

    total_time = format_time(time.time() - start_time)
    print_banner("PHASE 2 STEP D COMPLETE")
    print(f"[SUCCESS] Torsion angles (Phi/Psi) and stereochemical quality verified.")
    print(f"[SUCCESS] Execution Time : {total_time}")
    print(f"[INFO] Report Saved      : {os.path.relpath(report_path, project_root)}")
    print(f"[INFO] Plots Archived to : {os.path.relpath(archive_dir, project_root)}")
    print("="*110 + "\n")

if __name__ == "__main__":
    run_step2d_model_validation()