import os
import sys
import csv
import json
import time
import argparse
from datetime import datetime

# =============================================================================
# MINIMAL BOOTSTRAP -- locates the shared phase2_common module.
# See phase2_common.py for why this logic is centralized.
# =============================================================================
def _bootstrap_find_research_root(script_file):
    current = os.path.dirname(os.path.abspath(script_file))
    while os.path.basename(current) != "Research":
        parent = os.path.dirname(current)
        if parent == current:
            print(f"\n[FATAL ERROR] Could not locate a 'Research' anchor folder above: {script_file}")
            sys.exit(1)
        current = parent
    return current

_PROJECT_ROOT = _bootstrap_find_research_root(__file__)
_COMMON_DIR = os.path.join(_PROJECT_ROOT, "Phase 2", "_common")
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)

import phase2_common as common

# =============================================================================
# EXTERNALLY-SOURCED RESULTS -- MolProbity
#
# MolProbity has no simple pip/conda install and no public API. Rather than
# force a from-source build of the standalone Richardson Lab tool (risky on
# macOS) or a heavy full Phenix install just for one validation step, this
# uses the official free web server (molprobity.biochem.duke.edu) and reads
# back the values YOU record in manual_results.json -- the same pattern
# already used for Aggrescan3D, CamSol, and DeepSol in Step 2C.
#
# NOTE: Rampage is intentionally NOT used here. MolProbity's own summary
# report already includes Ramachandran favored/allowed/outlier percentages
# directly, so running a second, separate tool for the same numbers would
# be redundant.
# =============================================================================

MANUAL_RESULT_KEYS = {
    "molprobity_score": "MolProbity score (from the summary table)",
    "clashscore": "All-atom clashscore",
    "poor_rotamers_percent": "Poor rotamers (%)",
    "rama_favored_percent": "Ramachandran favored (%)",
    "rama_allowed_percent": "Ramachandran allowed (%)",
    "rama_outlier_percent": "Ramachandran outliers (%)",
}


def run_manual_prepare(structure_path, manual_results_path):
    common.print_banner("EXTERNAL RESULTS NEEDED -- MolProbity")
    print("Submit the structure below to the MolProbity web server, then fill")
    print("in the values in the template written below.")
    print("-" * 100)
    print("[MolProbity]  https://molprobity.biochem.duke.edu")
    print(f"              Upload: {structure_path}")
    print("              Run the full validation (default options are fine).")
    print("              Read the 6 values below off its summary table.")
    print("-" * 100)

    if os.path.isfile(manual_results_path):
        print(f"[INFO] {manual_results_path} already exists -- edit it directly, values are not overwritten.")
    else:
        template = {k: None for k in MANUAL_RESULT_KEYS}
        os.makedirs(os.path.dirname(manual_results_path), exist_ok=True)
        with open(manual_results_path, 'w') as f:
            json.dump(template, f, indent=2)
        print(f"[INFO] Template written to: {manual_results_path}")
        print("[INFO] Fill in each value after running MolProbity, then rerun this")
        print("[INFO] script normally (no --manual-prepare flag) to finish Step 2D.")
    print("=" * 100 + "\n")


def load_manual_results(manual_results_path):
    if not os.path.isfile(manual_results_path):
        return None
    with open(manual_results_path) as f:
        data = json.load(f)
    missing = [k for k, v in data.items() if v is None]
    if missing:
        return None
    return data

# =============================================================================
# CORE STEREOCHEMICAL VALIDATION DECISION ENGINE
# =============================================================================

def evaluate_model_quality(manual_results):
    """Evaluates stereochemical quality based on Phase II Step D rules."""

    result = {
        "MolProbity_Score": manual_results["molprobity_score"],
        "MolProbity_Clashscore": manual_results["clashscore"],
        "MolProbity_Poor_Rotamers_Pct": manual_results["poor_rotamers_percent"],
        "Rama_Favored": manual_results["rama_favored_percent"],
        "Rama_Allowed": manual_results["rama_allowed_percent"],
        "Rama_Outlier": manual_results["rama_outlier_percent"],
        "MP_Quality": "UNKNOWN",
        "Rama_Label": "UNKNOWN",
        "Rama_Meets_Target": False,
        "Overall_Status": "REVIEW",
    }

    # 1. MolProbity score bands, per methodology: 0.5-1.5 excellent,
    #    <2 acceptable, >3 poor.
    #    NOTE: the methodology text does not define the 2.0-3.0 range.
    #    Rather than silently folding that gap into either ACCEPTABLE or
    #    POOR, it is surfaced explicitly as MARGINAL so a real candidate
    #    landing there gets flagged for manual review instead of an
    #    unstated assumption deciding its fate.
    score = result["MolProbity_Score"]
    if 0.5 <= score <= 1.5:
        result["MP_Quality"] = "EXCELLENT"
    elif score < 2.0:
        result["MP_Quality"] = "ACCEPTABLE"
    elif score <= 3.0:
        result["MP_Quality"] = "MARGINAL (2.0-3.0 undefined by methodology)"
    else:
        result["MP_Quality"] = "POOR"

    # 2. Ramachandran: methodology's stated preferred target is 95-98%
    #    favored residues.
    favored = result["Rama_Favored"]
    result["Rama_Meets_Target"] = favored >= 95.0
    if 95.0 <= favored <= 98.0:
        result["Rama_Label"] = "MEETS TARGET (95-98%)"
    elif favored > 98.0:
        result["Rama_Label"] = "EXCEEDS TARGET (>98%)"
    else:
        result["Rama_Label"] = "BELOW TARGET (<95%)"

    # 3. Final verdict.
    #    NOTE: methodology grades MolProbity and Ramachandran independently
    #    but does not state how to combine them into one verdict. This
    #    requires BOTH to be satisfactory to call VALIDATED, treats a
    #    MARGINAL MolProbity score as REVIEW (not an automatic fail)
    #    rather than silently passing or failing it, and otherwise FAILS.
    if result["MP_Quality"] in ("EXCELLENT", "ACCEPTABLE") and result["Rama_Meets_Target"]:
        result["Overall_Status"] = "VALIDATED"
    elif "MARGINAL" in result["MP_Quality"]:
        result["Overall_Status"] = "REVIEW"
    else:
        result["Overall_Status"] = "FAILED - REQUIRES REFINEMENT"

    return result


def run_step2d_model_validation():
    start_time = time.time()
    project_root = _PROJECT_ROOT

    input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    # Validate the SAME minimized structure Step 2C analyzed -- not the raw
    # AlphaFold Server output from Step 2B, which is an earlier, unrepaired
    # version of the model that the rest of the pipeline no longer uses.
    stepc_archive_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC", "Supplementary_Archive")

    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepD")
    os.makedirs(output_base, exist_ok=True)

    common.print_banner("PHASE 2 STEP D: STEREOCHEMICAL MODEL VALIDATION")
    print(f"[INFO] Resolved Project Root : {project_root}")
    print("[INFO] Methodology : MolProbity (web server)")
    print("[INFO] Targets     : MolProbity 0.5-1.5 excellent / <2 acceptable / >3 poor | Favored 95-98%")
    print("-" * 110)

    winner_row, _ = common.get_winner_from_filtered_csv(input_csv_dir)
    if winner_row is None:
        return
    winner_name = winner_row["Variant"]
    safe_name = common.sanitize_variant_name(winner_name)

    structure_path = os.path.join(stepc_archive_dir, f"{safe_name}_minimized.pdb")
    if not os.path.isfile(structure_path):
        print(f"[ERROR] Structure file not found at {structure_path}.")
        print("[ERROR] Run Step 2C first so the minimized structure exists.")
        return

    print(f"[INFO] Validating Model : {os.path.relpath(structure_path, project_root)}")
    print("-" * 110)

    manual_results_path = os.path.join(output_base, f"{safe_name}_manual_results.json")
    manual_results = load_manual_results(manual_results_path)
    if manual_results is None:
        run_manual_prepare(structure_path, manual_results_path)
        return

    results = evaluate_model_quality(manual_results)

    print(f"{'VALIDATION METRIC':<30} | {'VALUE':<20} | {'ASSESSMENT'}")
    print("-" * 110)
    print(f"{'MolProbity Score':<30} | {results['MolProbity_Score']:<20.2f} | {results['MP_Quality']}")
    print(f"{'MolProbity Clashscore':<30} | {results['MolProbity_Clashscore']:<20.2f} | (contextual)")
    print(f"{'MolProbity Poor Rotamers %':<30} | {results['MolProbity_Poor_Rotamers_Pct']:<19.2f}% | (contextual)")
    print(f"{'Rama: Favored Region':<30} | {results['Rama_Favored']:<19.2f}% | {results['Rama_Label']}")
    print(f"{'Rama: Allowed Region':<30} | {results['Rama_Allowed']:<19.2f}% | (contextual)")
    print(f"{'Rama: Outlier Region':<30} | {results['Rama_Outlier']:<19.2f}% | (contextual)")
    print("-" * 110)
    print(f"FINAL DECISION : [{results['Overall_Status']}]")
    print("-" * 110)

    ts = datetime.now().strftime("%Y%m%d_%H%M")
    report_path = os.path.join(output_base, f"Step2D_Validation_Report_{ts}.csv")

    csv_data = {
        "Variant": winner_name,
        "Structure_File": os.path.relpath(structure_path, output_base),
        **results,
        "Manual_Results_File": os.path.relpath(manual_results_path, output_base),
    }

    with open(report_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=csv_data.keys())
        w.writeheader()
        w.writerow(csv_data)

    total_time = common.format_time(time.time() - start_time)
    common.print_banner("PHASE 2 STEP D COMPLETE")
    print("[SUCCESS] Stereochemical quality verified using real MolProbity results.")
    print(f"[SUCCESS] Execution Time : {total_time}")
    print(f"[INFO] Report Saved      : {os.path.relpath(report_path, project_root)}")
    print("=" * 110 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2D: Stereochemical Model Validation")
    parser.add_argument("--manual-prepare", action="store_true",
                         help="Print MolProbity submission instructions and write the results template early")
    args = parser.parse_args()

    if args.manual_prepare:
        project_root = _PROJECT_ROOT
        input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
        stepc_archive_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC", "Supplementary_Archive")
        output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepD")
        winner_row, _ = common.get_winner_from_filtered_csv(input_csv_dir)
        if winner_row:
            safe_name = common.sanitize_variant_name(winner_row["Variant"])
            structure_path = os.path.join(stepc_archive_dir, f"{safe_name}_minimized.pdb")
            manual_results_path = os.path.join(output_base, f"{safe_name}_manual_results.json")
            run_manual_prepare(structure_path, manual_results_path)
    else:
        run_step2d_model_validation()
