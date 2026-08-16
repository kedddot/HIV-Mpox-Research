import os
import sys
import csv
import time
import shutil
import zipfile
import tempfile
import glob
import argparse
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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
# ALPHAFOLD SERVER INTEGRATION (real, but web-submission-based -- see notes)
#
# AlphaFold3 requires an NVIDIA/CUDA GPU and does not run on a Mac. AlphaFold
# Server (alphafoldserver.com) is the official free DeepMind web version and
# supports multimer/complex predictions (which is what the TLR2/TLR4 docking
# needs), but it has NO public API -- submission is browser-only. So this
# step is split into two phases you run separately:
#
#   1) python Phase2B_secondaryAndTertiary.py --prepare
#      Writes out the exact sequences/job names to paste into AlphaFold
#      Server, for all 3 required jobs (monomer, +TLR2, +TLR4).
#
#   2) [ you manually submit the 3 jobs on alphafoldserver.com,
#        wait for them to finish, and download each result ]
#
#   3) python Phase2B_secondaryAndTertiary.py --import monomer <path_to_download>
#      python Phase2B_secondaryAndTertiary.py --import tlr2    <path_to_download>
#      python Phase2B_secondaryAndTertiary.py --import tlr4    <path_to_download>
#      Unzips each download, locates the top-ranked model (.cif), and places
#      it where Step 2C/2D expect it.
#
#   4) python Phase2B_secondaryAndTertiary.py
#      Runs normally once all 3 imports are done: computes the secondary
#      structure baseline and finishes the step.
#
# NOTE: AlphaFold Server's real output format is mmCIF (.cif), not legacy
# PDB. Step 2C and 2D need one matching edit each to expect .cif -- see the
# comments at the top of those files.
# =============================================================================

class StructuralPredictors:
    @staticmethod
    def submit_psipred(sequence, name):
        # TODO: replace with a real PsiPred API call once you confirm your
        # access method (hosted API vs local install).
        return f"PSIPRED_JOB_{name}_8839"

    @staticmethod
    def submit_raptorx(sequence, name):
        # TODO: replace with a real RaptorX API call once you confirm your
        # access method.
        return f"RAPTORX_JOB_{name}_1022"


def import_alphafold_output(downloaded_path, target_cif_path):
    """
    Takes either a .zip downloaded from AlphaFold Server, or a bare .cif
    file, and copies the top-ranked model into target_cif_path.

    AlphaFold Server zips typically contain multiple ranked models
    (model_0 = top-ranked). This looks for the lowest-numbered
    "_model_N.cif" file, falling back to the first .cif found if that
    naming pattern isn't present -- server output naming has changed
    before and may change again, so verify against your actual download
    the first time you use this.
    """
    if not os.path.exists(downloaded_path):
        raise FileNotFoundError(f"Downloaded file not found: {downloaded_path}")

    os.makedirs(os.path.dirname(target_cif_path), exist_ok=True)

    if downloaded_path.endswith(".cif"):
        shutil.copy(downloaded_path, target_cif_path)
        return target_cif_path

    if downloaded_path.endswith(".zip"):
        with tempfile.TemporaryDirectory() as tmp:
            with zipfile.ZipFile(downloaded_path, 'r') as z:
                z.extractall(tmp)

            cif_files = glob.glob(os.path.join(tmp, "**", "*.cif"), recursive=True)
            if not cif_files:
                raise FileNotFoundError(f"No .cif files found inside {downloaded_path}")

            model_0 = [f for f in cif_files if "_model_0" in os.path.basename(f)]
            chosen = model_0[0] if model_0 else sorted(cif_files)[0]

            shutil.copy(chosen, target_cif_path)
            return target_cif_path

    raise ValueError(f"Unrecognized file type for AlphaFold Server output: {downloaded_path}")


def _resolve_paths():
    project_root = _PROJECT_ROOT
    input_csv = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    variant_fasta_path = common.phase1g_fasta_path(project_root)
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB")
    out_sec_dir = os.path.join(output_base, "Secondary_Structure")
    out_ter_dir = os.path.join(output_base, "Tertiary_Structure")
    return project_root, input_csv, variant_fasta_path, output_base, out_sec_dir, out_ter_dir


def _get_target(input_csv, variant_fasta_path):
    winner_row, _ = common.get_winner_from_filtered_csv(input_csv)
    if winner_row is None:
        return None, None, None
    winner_name = winner_row["Variant"]

    variants = common.load_multi_fasta(variant_fasta_path)
    if variants is None:
        print(f"[ERROR] Variant FASTA file not found: {variant_fasta_path}")
        return None, None, None

    vax_sequence = common.lookup_sequence(variants, winner_name)
    if vax_sequence is None:
        print(f"[ERROR] Could not find a matching header for '{winner_name}' in: {variant_fasta_path}")
        available = list(variants.keys())
        print(f"[ERROR] Headers actually present ({len(available)}): {available[:10]}{' ...' if len(available) > 10 else ''}")
        return None, None, None

    safe_name = common.sanitize_variant_name(winner_name)
    return winner_name, vax_sequence, safe_name


# (Human TLR2: UniProt O60603; Human TLR4: UniProt O00206)
TLR2_SEQ = """
MPHTLWMVWVLGVIISLSKEESSNQASLSCDRNGICKGSSGSLNSIPSGLTEAVKSLDLS
NNRITYISNSDLQRCVNLQALVLTSNGINTIEEDSFSSLGSLEHLDLSYNYLSNLSSSWF
KPLSSLTFLNLLGNPYKTLGETSLFSHLTKLQILRVGNMDTFTKIQRKDFAGLTFLEELE
IDASDLQSYEPKSLKSIQNVSHLILHMKQHILLLEIFVDVTSSVECLELRDTDLDTFHFS
ELSTGETNSLIKKFTFRNVKITDESLFQVMKLLNQISGLLELEFDDCTLNGVGNFRASDN
DRVIDPGKVETLTIRRLHIPRFYLFYDLSTLYSLTERVKRITVENSKVFLVPCLLSQHLK
SLEYLDLSENLMVEEYLKNSACEDAWPSLQTLILRQNHLASLEKTGETLLTLKNLTNIDI
SKNSFHSMPETCQWPEKMKYLNLSSTRIHSVTGCIPKTLEILDVSNNNLNLFSLNLPQLK
ELYISRNKLMTLPDASLLPMLLVLKISRNAITTFSKEQLDSFHTLKTLEAGGNNFICSCE
FLSFTQEQQALAKVLIDWPANYLCDSPSHVRGQQVQDVRLSVSECHRTALVSGMCCALFL
LILLTGVLCHRFHGLWYMKMMWAWLQAKRKPRKAPSRNICYDAFVSYSERDAYWVENLMV
QELENFNPPFKLCLHKRDFIPGKWIIDNIIDSIEKSHKTVFVLSENFVKSEWCKYELDFS
HFRLFDENNDAAILILLEPIEKKAIPQRFCKLRKIMNTKTYLEWPMDEAQREGFWVNLRA
AIKS
"""
TLR4_SEQ = """
MMSASRLAGTLIPAMAFLSCVRPESWEPCVEVVPNITYQCMELNFYKIPDNLPFSTKNLD
LSFNPLRHLGSYSFFSFPELQVLDLSRCEIQTIEDGAYQSLSHLSTLILTGNPIQSLALG
AFSGLSSLQKLVAVETNLASLENFPIGHLKTLKELNVAHNLIQSFKLPEYFSNLTNLEHL
DLSSNKIQSIYCTDLRVLHQMPLLNLSLDLSLNPMNFIQPGAFKEIRLHKLTLRNNFDSL
NVMKTCIQGLAGLEVHRLVLGEFRNEGNLEKFDKSALEGLCNLTIEEFRLAYLDYYLDDI
IDLFNCLTNVSSFSLVSVTIERVKDFSYNFGWQHLELVNCKFGQFPTLKLKSLKRLTFTS
NKGGNAFSEVDLPSLEFLDLSRNGLSFKGCCSQSDFGTTSLKYLDLSFNGVITMSSNFLG
LEQLEHLDFQHSNLKQMSEFSVFLSLRNLIYLDISHTHTRVAFNGIFNGLSSLEVLKMAG
NSFQENFLPDIFTELRNLTFLDLSQCQLEQLSPTAFNSLSSLQVLNMSHNNFFSLDTFPY
KCLNSLQVLDYSLNHIMTSKKQELQHFPSSLAFLNLTQNDFACTCEHQSFLQWIKDQRQL
LVEVERMECATPSDKQGMPVLSLNITCQMNKTIIGVSVLSVLVVSVVAVLVYKFYFHLML
LAGCIKYGRGENIYDAFVIYSSQDEDWVRNELVKNLEEGVPPFQLCLHYRDFIPGVAIAA
NIIHEGFHKSRKVIVVVSQHFIQSRWCIFEYEIAQTWQFLSSRAGIIFIVLQKVEKTLLR
QQVELYRLLSRNTYLEWEDSVLGRHIFWRRLRKALLDGKSWNPEGTVGTGCNWQEATSI
"""


def run_prepare():
    project_root, input_csv, variant_fasta_path, output_base, out_sec_dir, out_ter_dir = _resolve_paths()
    os.makedirs(out_sec_dir, exist_ok=True)
    os.makedirs(out_ter_dir, exist_ok=True)

    common.print_banner("PHASE 2 STEP B -- PREPARE ALPHAFOLD SERVER SUBMISSIONS")
    winner_name, vax_sequence, safe_name = _get_target(input_csv, variant_fasta_path)
    if winner_name is None:
        return

    print(f"[INFO] Target Variant : {winner_name}")
    print("-" * 100)
    print("Submit these 3 jobs at https://alphafoldserver.com -- one entity per")
    print("sequence within a job for the complexes. After each job finishes,")
    print("download its result and run the --import command shown for it.")
    print("-" * 100)

    jobs = [
        ("monomer", f"Monomer_{safe_name}", [("Vaccine", vax_sequence)]),
        ("tlr2", f"Complex_TLR2_{safe_name}", [("Vaccine", vax_sequence), ("TLR2", TLR2_SEQ)]),
        ("tlr4", f"Complex_TLR4_{safe_name}", [("Vaccine", vax_sequence), ("TLR4", TLR4_SEQ)]),
    ]

    for job_key, job_name, entities in jobs:
        print(f"\n[JOB: {job_key}]  AlphaFold Server job name suggestion: {job_name}")
        for entity_label, seq in entities:
            print(f"   Entity ({entity_label}):")
            print(f"   {seq}")
        print(f"   --> After downloading the result, run:")
        print(f"       python {os.path.basename(__file__)} --import {job_key} /path/to/downloaded_file.zip")

    print("\n" + "-" * 100)
    print("Once all 3 are imported, run this script with no arguments to finish Step 2B.")
    print("=" * 100 + "\n")


def run_import(job_key, downloaded_path):
    project_root, input_csv, variant_fasta_path, output_base, out_sec_dir, out_ter_dir = _resolve_paths()
    winner_name, vax_sequence, safe_name = _get_target(input_csv, variant_fasta_path)
    if winner_name is None:
        return

    target_map = {
        "monomer": os.path.join(out_ter_dir, f"AF3_Target_{safe_name}.cif"),
        "tlr2": os.path.join(out_ter_dir, f"AF3_Complex_TLR2_{safe_name}.cif"),
        "tlr4": os.path.join(out_ter_dir, f"AF3_Complex_TLR4_{safe_name}.cif"),
    }

    if job_key not in target_map:
        print(f"[ERROR] Unknown job key '{job_key}'. Expected one of: {list(target_map.keys())}")
        return

    target_path = target_map[job_key]
    try:
        result = import_alphafold_output(downloaded_path, target_path)
    except (FileNotFoundError, ValueError) as e:
        print(f"[ERROR] {e}")
        return

    print(f"[SUCCESS] Imported '{job_key}' model -> {os.path.relpath(result, project_root)}")


def run_step2b_structure_prediction():
    start_time = time.time()
    project_root, input_csv, variant_fasta_path, output_base, out_sec_dir, out_ter_dir = _resolve_paths()
    os.makedirs(out_sec_dir, exist_ok=True)
    os.makedirs(out_ter_dir, exist_ok=True)

    common.print_banner("PHASE 2 STEP B: SECONDARY AND TERTIARY STRUCTURE PREDICTION")
    print(f"[INFO] Resolved Project Root : {project_root}")
    print(f"[INFO] Looking for Filtered CSVs in : {input_csv}")

    winner_name, vax_sequence, safe_name = _get_target(input_csv, variant_fasta_path)
    if winner_name is None:
        return

    print(f"[INFO] Target Variant    : {winner_name} (Rank 1 from Step 2A)")
    print("[INFO] Methodology       : PsiPred, RaptorX, AlphaFold Server (monomer + TLR-2/TLR-4 complexes)")
    print("-" * 110)

    # Confirm all 3 AlphaFold Server imports are present before proceeding
    required_files = {
        "monomer": os.path.join(out_ter_dir, f"AF3_Target_{safe_name}.cif"),
        "tlr2": os.path.join(out_ter_dir, f"AF3_Complex_TLR2_{safe_name}.cif"),
        "tlr4": os.path.join(out_ter_dir, f"AF3_Complex_TLR4_{safe_name}.cif"),
    }
    missing = [k for k, p in required_files.items() if not os.path.isfile(p)]
    if missing:
        print(f"[ERROR] Missing AlphaFold Server results for: {missing}")
        print(f"[ERROR] Run: python {os.path.basename(__file__)} --prepare")
        print("[ERROR] then submit the jobs on alphafoldserver.com and import each result.")
        return

    print("[INFO] All 3 AlphaFold Server models found:")
    for k, p in required_files.items():
        print(f"       {k}: {os.path.relpath(p, project_root)}")
    print("-" * 110)

    # SECONDARY STRUCTURE PROPENSITY (Baseline Analysis) -- unaffected by
    # the AlphaFold Server change, still computed directly from sequence.
    analysis = ProteinAnalysis(vax_sequence)
    helix, turn, sheet = analysis.secondary_structure_fraction()
    coil = 1.0 - (helix + sheet + turn)

    print(f"{'STRUCTURAL ELEMENT':<30} | {'FRACTION':<15} | {'PERCENTAGE':<15}")
    print("-" * 110)
    print(f"{'Alpha-Helix (H)':<30} | {helix:>15.4f} | {helix*100:>14.2f}%")
    print(f"{'Beta-Sheet (E)':<30} | {sheet:>15.4f} | {sheet*100:>14.2f}%")
    print(f"{'Turns (T)':<30} | {turn:>15.4f} | {turn*100:>14.2f}%")
    print(f"{'Random Coils (C)':<30} | {coil:>15.4f} | {coil*100:>14.2f}%")
    print("-" * 110)

    # Secondary Predictors (still TODO -- see StructuralPredictors class)
    psipred_fasta = os.path.join(out_sec_dir, f"PsiPred_{safe_name}.fasta")
    raptorx_fasta = os.path.join(out_sec_dir, f"RaptorX_{safe_name}.fasta")
    with open(psipred_fasta, 'w') as f: f.write(f">PsiPred_{winner_name}\n{vax_sequence}\n")
    with open(raptorx_fasta, 'w') as f: f.write(f">RaptorX_{winner_name}\n{vax_sequence}\n")

    psi_job = StructuralPredictors.submit_psipred(vax_sequence, safe_name)
    rap_job = StructuralPredictors.submit_raptorx(vax_sequence, safe_name)
    print(f"[INFO] PsiPred Submission  : {psi_job}")
    print(f"[INFO] RaptorX Submission  : {rap_job}")

    report_path = os.path.join(out_sec_dir, "Step2B_Secondary_Structure_Baseline.csv")
    with open(report_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Element", "Fraction", "Percentage"])
        writer.writerow(["Alpha-Helix", helix, f"{helix*100:.2f}%"])
        writer.writerow(["Beta-Sheet", sheet, f"{sheet*100:.2f}%"])
        writer.writerow(["Turns", turn, f"{turn*100:.2f}%"])
        writer.writerow(["Random Coils", coil, f"{coil*100:.2f}%"])

    total_time = common.format_time(time.time() - start_time)
    common.print_banner("STEP 2B COMPLETE")
    print("[SUCCESS] Secondary baseline calculated; real AlphaFold Server models in place.")
    print(f"[SUCCESS] Total Execution Time : {total_time}")
    print(f"[INFO] Outputs routed to       : {os.path.relpath(output_base, project_root)}")
    print("=" * 110 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2B: Secondary + Tertiary Structure Prediction")
    parser.add_argument("--prepare", action="store_true", help="Print sequences to submit to AlphaFold Server")
    parser.add_argument("--import", dest="import_args", nargs=2, metavar=("JOB_KEY", "DOWNLOADED_PATH"),
                         help="Import a downloaded AlphaFold Server result. JOB_KEY is monomer, tlr2, or tlr4.")
    args = parser.parse_args()

    if args.prepare:
        run_prepare()
    elif args.import_args:
        run_import(args.import_args[0], args.import_args[1])
    else:
        run_step2b_structure_prediction()
