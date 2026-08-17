import os
import sys
import csv
import json
import time
import shutil
import subprocess
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
# CONFIGURE -- fill these in for your machine before running
# =============================================================================
FOLDX_BINARY = "/Users/nek/Desktop/School/Research/FoldX5/foldx5.1_20261231"
APBS_BINARY = "/opt/miniconda3/envs/apbs_env/bin/apbs"

# =============================================================================
# REAL TOOLS -- repair, protonation, SASA, FoldX, APBS
# =============================================================================

def repair_structure_pdbfixer(input_path, output_path):
    """
    Structure repair (methodology prep step 1, part A) via OpenMM PDBFixer.
    Adds missing atoms/residues, removes heterogens.

    NOTE: this step alone does NOT relax/minimize the structure -- it only
    patches in missing atoms. See minimize_structure_openmm() below for the
    actual "minimal relaxation" half of prep step 1; both must run in
    sequence, or FoldX (and to a lesser extent SASA) can report wildly
    unrealistic values off of un-relaxed steric clashes.

    Requires: pip install pdbfixer openmm   (or via conda-forge if pip
    gives you trouble with the OpenMM C-extension dependency)
    """
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=input_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.removeHeterogens(keepWater=False)

    with open(output_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return output_path


def minimize_structure_openmm(input_pdb_path, output_pdb_path, max_iterations=500):
    """
    Real energy minimization (methodology prep step 1, part B) via OpenMM,
    using the AMBER14 forcefield with GBn2 implicit solvent. This resolves
    local steric clashes left over from PDBFixer's atom-repair alone --
    AlphaFold-predicted multi-domain constructs joined by short linkers are
    especially prone to minor clashes right at the junctions, which
    downstream tools (particularly FoldX) are very sensitive to: an
    un-minimized structure can report wildly unrealistic large-positive
    ddG values that reflect clash artifacts, not genuine instability.

    Requires: pip install openmm (already installed alongside pdbfixer)
    """
    from openmm.app import PDBFile, ForceField, Modeller, Simulation, HBonds, NoCutoff
    from openmm import LangevinMiddleIntegrator
    from openmm.unit import kelvin, picosecond, picoseconds

    pdb = PDBFile(input_pdb_path)
    forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml')

    modeller = Modeller(pdb.topology, pdb.positions)
    # PDBFixer's repair step only adds missing HEAVY atoms -- it never adds
    # hydrogens, so the forcefield below (which needs a fully complete
    # molecule to match its templates) can't build a System without them
    # yet. This is a rough, generic hydrogen placement purely so
    # minimization has something physically complete to work with -- it is
    # NOT the rigorous, pKa-aware pH 7.0 protonation your methodology
    # calls for; that happens later, more accurately, via PDB2PQR/PROPKA.
    modeller.addHydrogens(forcefield, pH=7.0)
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
    )
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    simulation.minimizeEnergy(maxIterations=max_iterations)
    positions = simulation.context.getState(getPositions=True).getPositions()

    # Strip hydrogens back out before writing the final file. They were
    # necessary for the forcefield's physics during minimization, but
    # every downstream consumer of this file -- our own FreeSASA call,
    # CamSol, Aggrescan3D, and PDB2PQR's own re-protonation -- expects a
    # heavy-atom-only structure, matching virtually any real PDB/
    # AlphaFold structure. Leaving hydrogens in caused external tools
    # (e.g. CamSol's own FreeSASA-based backend) to crash outright, and
    # likely silently skewed our own local SASA numbers too.
    minimized_modeller = Modeller(modeller.topology, positions)
    hydrogens = [a for a in minimized_modeller.topology.atoms() if a.element is not None and a.element.symbol == 'H']
    minimized_modeller.delete(hydrogens)

    with open(output_pdb_path, 'w') as f:
        PDBFile.writeFile(minimized_modeller.topology, minimized_modeller.positions, f)
    return output_pdb_path


def protonate_pdb2pqr(input_pdb_path, output_pqr_path, ph=7.0):
    """
    Real protonation at target pH (methodology prep step 2) via PDB2PQR.
    Requires: pip install pdb2pqr
    """
    cmd = [
        "pdb2pqr30",
        "--ff=AMBER",
        f"--with-ph={ph}",
        "--titration-state-method=propka",
        input_pdb_path,
        output_pqr_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"PDB2PQR failed (exit {result.returncode}):\n{result.stderr}")
    return output_pqr_path


HYDROPHOBIC_RESIDUES = {'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'MET', 'TRP', 'CYS'}
EXPOSURE_THRESHOLD_A2 = 15.0  # per-residue SASA above this counts as "exposed"


def analyze_sasa_freesasa(pdb_path):
    """
    Real SASA + hydrophobic surface fraction via FreeSASA.
    Requires: pip install freesasa

    NOTE on "exposed_patch_area": true spatial patch detection (like
    Aggrescan3D does) clusters residues by 3D proximity. This computes a
    simpler, sequence-adjacency-based approximation -- the longest run of
    consecutive, exposed, hydrophobic residues along the chain. That is a
    real, defensible metric, but it is NOT equivalent to a true 3D
    spatial-clustering patch. Flagged here so you don't mistake this for
    Aggrescan3D-grade analysis; Aggrescan3D itself is still handled
    manually (see run_step2c_manual_prepare / manual_results.json).
    """
    import freesasa

    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    residue_areas = result.residueAreas()

    total_sasa = result.totalArea()
    hydrophobic_sasa = 0.0
    max_patch_area = 0.0
    current_patch = 0.0

    for chain_id in residue_areas:
        for res_num in sorted(residue_areas[chain_id], key=lambda x: int(x)):
            area = residue_areas[chain_id][res_num]
            is_hydrophobic = area.residueType in HYDROPHOBIC_RESIDUES
            is_exposed = area.total > EXPOSURE_THRESHOLD_A2

            if is_hydrophobic:
                hydrophobic_sasa += area.total

            if is_hydrophobic and is_exposed:
                current_patch += area.total
                max_patch_area = max(max_patch_area, current_patch)
            else:
                current_patch = 0.0

    hydrophobic_fraction = hydrophobic_sasa / total_sasa if total_sasa > 0 else 0.0

    return {
        "hydrophobic_fraction": hydrophobic_fraction,
        "exposed_patch_area": max_patch_area,
        "total_sasa": total_sasa,
    }


def run_foldx(pdb_path, output_dir):
    """
    Real FoldX Stability command -> ddG.
    Requires: FoldX binary, free academic registration at foldxsuite.crg.eu
    Set FOLDX_BINARY at the top of this file to your install path.

    NOTE: FoldX's exact CLI flags and output filename have changed across
    versions -- verify the "_0_ST.fxout" output name matches your version
    the first time you run this.
    """
    if not os.path.isfile(FOLDX_BINARY):
        raise FileNotFoundError(
            f"FoldX binary not found at '{FOLDX_BINARY}'. Set FOLDX_BINARY "
            f"at the top of this file to your actual FoldX install path."
        )

    pdb_dir = os.path.dirname(pdb_path)
    pdb_name = os.path.basename(pdb_path)

    # FoldX reads rotabase.txt from its OWN working directory at runtime,
    # but we run it with cwd set to the PDB's folder (below) -- so copy
    # rotabase.txt alongside the PDB each time rather than relying on you
    # to remember to place it there manually. It must live next to the
    # FoldX binary itself (that's the standard FoldX distribution layout).
    rotabase_src = os.path.join(os.path.dirname(FOLDX_BINARY), "rotabase.txt")
    rotabase_dst = os.path.join(pdb_dir, "rotabase.txt")
    if os.path.isfile(rotabase_src) and not os.path.isfile(rotabase_dst):
        shutil.copy(rotabase_src, rotabase_dst)
    elif not os.path.isfile(rotabase_src):
        raise FileNotFoundError(
            f"rotabase.txt not found next to FOLDX_BINARY at '{rotabase_src}'. "
            f"FoldX needs this file to run -- it should have come with your "
            f"FoldX download; place it in the same folder as the binary."
        )

    cmd = [
        FOLDX_BINARY,
        "--command=Stability",
        f"--pdb={pdb_name}",
        f"--output-dir={output_dir}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=pdb_dir)
    if result.returncode != 0:
        raise RuntimeError(f"FoldX failed (exit {result.returncode}):\n{result.stderr}")

    stability_file = os.path.join(output_dir, f"{os.path.splitext(pdb_name)[0]}_0_ST.fxout")
    if not os.path.isfile(stability_file):
        raise FileNotFoundError(
            f"Expected FoldX output not found: {stability_file}\n"
            f"Your FoldX version may name its output differently -- check {output_dir} manually."
        )

    with open(stability_file) as f:
        line = f.readline()
    ddg = float(line.strip().split("\t")[1])
    return ddg


def run_apbs(pqr_path, output_dir):
    """
    Real APBS electrostatic surface map.
    Requires: conda install -c bioconda -c conda-forge apbs
    Set APBS_BINARY at the top of this file if 'apbs' isn't on your PATH.

    NOTE: grid dimensions below (dime/cglen/fglen) are generic defaults.
    For a production run, tune these to your structure's actual size --
    see APBS documentation for guidance on grid spacing vs. structure extent.
    """
    if shutil.which(APBS_BINARY) is None and not os.path.isfile(APBS_BINARY):
        raise FileNotFoundError(
            f"APBS binary not found ('{APBS_BINARY}'). Set APBS_BINARY at "
            f"the top of this file, or ensure 'apbs' is on your PATH."
        )

    apbs_input_path = os.path.join(output_dir, "apbs_input.in")
    dx_output_prefix = os.path.join(output_dir, "electrostatic_surface")

    apbs_input = f"""read
    mol pqr {pqr_path}
end
elec
    mg-auto
    dime 97 97 97
    cglen 100 100 100
    fglen 80 80 80
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.54
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx {dx_output_prefix}
end
quit
"""
    with open(apbs_input_path, 'w') as f:
        f.write(apbs_input)

    cmd = [APBS_BINARY, apbs_input_path]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)
    if result.returncode != 0:
        raise RuntimeError(f"APBS failed (exit {result.returncode}):\n{result.stderr}")

    return f"{dx_output_prefix}.dx"


# =============================================================================
# EXTERNALLY-SOURCED RESULTS -- Aggrescan3D, CamSol (structural mode), DeepSol
#
# None of these have a free, one-command-automatable path from this script.
# Rather than fake them with more mocks, this reads real values YOU obtain
# separately and record in manual_results.json.
#
#   - Aggrescan3D: run via the LOCAL standalone install, in your own
#     terminal, separately from this script (not the aggrescan3d.pl web
#     server). Record its results the same way as the other two below.
#   - CamSol (structural mode) and DeepSol: web-submission tools, no
#     local install used.
# =============================================================================

MANUAL_RESULT_KEYS = {
    "aggrescan3d_max_patch_length": "Aggrescan3D max contiguous aggregation-prone patch length (aa)",
    "aggrescan3d_in_high_plddt_region": "Whether that patch falls in a pLDDT>=70 region (true/false)",
    "camsol_structural_max_patch_length": "CamSol (structural mode) max aggregation patch length (aa)",
    "camsol_structural_in_high_plddt_region": "Whether that patch falls in a pLDDT>=70 region (true/false)",
    # NOTE: SOLpro's server is no longer reachable, so this pipeline now uses
    # DeepSol S2 (Khurana et al. 2018) as the sequence-based solubility
    # predictor instead. DeepSol also outputs a genuine 0.0-1.0 probability
    # (unlike e.g. CamSol Intrinsic's unbounded score), so the decision rule
    # below (`< 0.50`) is still valid unchanged -- only the source tool and
    # field name changed.
    "deepsol_probability": "DeepSol S2 predicted solubility probability (0.0-1.0)",
}


def run_manual_prepare(sequence, repaired_pdb_path, manual_results_path):
    common.print_banner("EXTERNAL RESULTS NEEDED -- Aggrescan3D, CamSol, DeepSol")
    print("Aggrescan3D: run the LOCAL standalone install in your own terminal")
    print("(not the aggrescan3d.pl web server) on the repaired PDB below.")
    print("CamSol and DeepSol: submit manually through each site.")
    print("Then fill in the values in the template written below.")
    print("-" * 100)

    # Write the FASTA ourselves from the same in-pipeline sequence used
    # everywhere else in Step C, rather than relying on a hand-copied
    # sequence string -- this is the exact sequence loaded from the
    # Phase 1G FASTA earlier in run_step2c_solubility_analysis().
    archive_dir = os.path.dirname(manual_results_path)
    fasta_path = os.path.join(archive_dir, f"{os.path.basename(repaired_pdb_path).replace('_repaired.pdb', '')}.fasta")
    variant_id = os.path.basename(fasta_path).replace(".fasta", "")
    os.makedirs(archive_dir, exist_ok=True)
    with open(fasta_path, 'w') as f:
        f.write(f">{variant_id}\n{sequence}\n")
    print(f"[SUCCESS] FASTA written: {fasta_path}")
    print("-" * 100)

    print(f"[Aggrescan3D]  Local standalone install -- run on: {repaired_pdb_path}")
    print(f"[CamSol]       https://www-cohsoftware.ch.cam.ac.uk/index.php/camsolstructural -- upload: {repaired_pdb_path}")
    print(f"[DeepSol]      https://machinelearning-protein.qcri.org  -- upload: {fasta_path}")
    print(f"               (DeepSol replaces SOLpro -- SOLpro's server is no longer reachable.")
    print(f"               On the DeepSol form, set the model parameter to 2 -- DeepSol S2 was")
    print(f"               the best-performing of the three published architectures in the")
    print(f"               original paper, Khurana et al. 2018. NOTE: DeepSol's web server")
    print(f"               requires a free QCAI account before you can submit a job -- sign")
    print(f"               up first if you haven't already. If this URL doesn't work, try")
    print(f"               qcai.qcri.org instead -- QCRI has been consolidating tools there.)")
    print("-" * 100)

    if os.path.isfile(manual_results_path):
        print(f"[INFO] {manual_results_path} already exists -- edit it directly, values are not overwritten.")
    else:
        template = {k: None for k in MANUAL_RESULT_KEYS}
        os.makedirs(os.path.dirname(manual_results_path), exist_ok=True)
        with open(manual_results_path, 'w') as f:
            json.dump(template, f, indent=2)
        print(f"[INFO] Template written to: {manual_results_path}")
        print("[INFO] Fill in each value after running the tools above, then rerun this")
        print("[INFO] script normally (no --manual-prepare flag) to finish Step 2C.")
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
# SOLUBILITY & STRUCTURAL INTEGRITY DECISION ENGINE
# =============================================================================

def evaluate_structural_solubility(sequence, repaired_pdb_path, sasa_result, ddg, manual_results):
    deepsol_prob = manual_results["deepsol_probability"]
    agg_patch_len = manual_results["aggrescan3d_max_patch_length"]
    agg_high_plddt = manual_results["aggrescan3d_in_high_plddt_region"]
    camsol_patch_len = manual_results["camsol_structural_max_patch_length"]
    camsol_high_plddt = manual_results["camsol_structural_in_high_plddt_region"]

    result = {
        "DeepSol_Prob": deepsol_prob,
        "Agg_Patch_Len": agg_patch_len,
        "Agg_High_pLDDT": agg_high_plddt,
        "CamSol_Patch_Len": camsol_patch_len,
        "CamSol_High_pLDDT": camsol_high_plddt,
        "Hydrophobic_Fraction": sasa_result["hydrophobic_fraction"],
        "Exposed_Area": sasa_result["exposed_patch_area"],
        "FoldX_ddG": ddg,
        "Status": "PASS",
        "Reason_Code": [],
    }

    is_insoluble = deepsol_prob < 0.50
    has_confident_agg_patch = agg_patch_len > 8 and agg_high_plddt

    if is_insoluble and has_confident_agg_patch:
        result["Status"] = "REJECT"
        result["Reason_Code"].append("DeepSol <0.5 & Aggrescan3D Confident Patch >8aa")

    if result["Status"] != "REJECT":
        if sasa_result["hydrophobic_fraction"] > 0.25:
            result["Status"] = "REVIEW"
            result["Reason_Code"].append("Hydrophobic Fraction > 0.25")
        if sasa_result["exposed_patch_area"] > 250.0:
            result["Status"] = "REVIEW"
            result["Reason_Code"].append("Exposed Hydrophobic Area > 250 \u00c5\u00b2")
        if ddg > 1.5:
            result["Status"] = "REVIEW"
            result["Reason_Code"].append("FoldX \u0394\u0394G > +1.5 kcal/mol")
        if camsol_patch_len > 8 and camsol_high_plddt:
            result["Status"] = "REVIEW"
            result["Reason_Code"].append("CamSol (structural) Aggregation Patch > 8aa")

    if not result["Reason_Code"]:
        result["Reason_Code"].append("All structural parameters optimal")

    result["Reason_Code"] = " | ".join(result["Reason_Code"])
    return result


def run_step2c_solubility_analysis():
    start_time = time.time()
    project_root = _PROJECT_ROOT

    input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    variant_fasta_path = common.phase1g_fasta_path(project_root)
    pdb_input_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB", "Tertiary_Structure")

    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC")
    archive_dir = os.path.join(output_base, "Supplementary_Archive")
    os.makedirs(output_base, exist_ok=True)
    os.makedirs(archive_dir, exist_ok=True)

    common.print_banner("PHASE 2 STEP C: 3D SOLUBILITY & STRUCTURAL INTEGRITY")
    print(f"[INFO] Resolved Project Root : {project_root}")

    winner_row, _ = common.get_winner_from_filtered_csv(input_csv_dir)
    if winner_row is None:
        return
    winner_name = winner_row["Variant"]

    variants = common.load_multi_fasta(variant_fasta_path)
    if variants is None:
        print(f"[ERROR] Variant FASTA file not found: {variant_fasta_path}")
        return
    sequence = common.lookup_sequence(variants, winner_name)
    if sequence is None:
        print(f"[ERROR] Could not find a matching header for '{winner_name}' in: {variant_fasta_path}")
        return

    safe_name = common.sanitize_variant_name(winner_name)
    input_cif_path = os.path.join(pdb_input_dir, f"AF3_Target_{safe_name}.cif")

    if not os.path.isfile(input_cif_path):
        print(f"[ERROR] Structure file not found at {input_cif_path}.")
        print("[ERROR] Run Step 2B first (including the AlphaFold Server import) so this model exists.")
        return

    print(f"[INFO] Target Variant : {winner_name}")
    print(f"[INFO] Input Structure: {input_cif_path}")
    print("-" * 110)

    # 1a. Repair (real, PDBFixer) -- adds missing atoms/residues only
    repaired_pdb_path = os.path.join(archive_dir, f"{safe_name}_repaired.pdb")
    try:
        repair_structure_pdbfixer(input_cif_path, repaired_pdb_path)
        print(f"[SUCCESS] Repaired structure   : {os.path.relpath(repaired_pdb_path, output_base)}")
    except ImportError:
        print("[ERROR] PDBFixer/OpenMM not installed. Run: pip install pdbfixer openmm")
        return
    except Exception as e:
        print(f"[ERROR] PDBFixer repair failed: {e}")
        return

    # 1b. Minimize (real, OpenMM) -- resolves steric clashes left over from
    # repair alone. All downstream steps use this minimized structure, not
    # the raw-repaired one, since un-minimized clashes can otherwise blow
    # up FoldX ddG into unrealistic large-positive values.
    minimized_pdb_path = os.path.join(archive_dir, f"{safe_name}_minimized.pdb")
    try:
        minimize_structure_openmm(repaired_pdb_path, minimized_pdb_path)
        print(f"[SUCCESS] Minimized structure   : {os.path.relpath(minimized_pdb_path, output_base)}")
    except ImportError:
        print("[ERROR] OpenMM not installed. Run: pip install openmm")
        return
    except Exception as e:
        print(f"[ERROR] OpenMM minimization failed: {e}")
        return
    repaired_pdb_path = minimized_pdb_path  # everything downstream uses the minimized structure

    # 2. Protonate at pH 7.0 (real, PDB2PQR)
    protonated_pqr_path = os.path.join(archive_dir, f"{safe_name}_pH7.0.pqr")
    try:
        protonate_pdb2pqr(repaired_pdb_path, protonated_pqr_path, ph=7.0)
        print(f"[SUCCESS] Protonated (pH 7.0): {os.path.relpath(protonated_pqr_path, output_base)}")
    except FileNotFoundError:
        print("[ERROR] pdb2pqr30 not found on PATH. Run: pip install pdb2pqr")
        return
    except Exception as e:
        print(f"[ERROR] PDB2PQR protonation failed: {e}")
        return

    # 3. SASA + hydrophobic fraction (real, FreeSASA)
    try:
        sasa_result = analyze_sasa_freesasa(repaired_pdb_path)
        print(f"[SUCCESS] SASA analysis complete (hydrophobic fraction: {sasa_result['hydrophobic_fraction']:.3f})")
    except ImportError:
        print("[ERROR] freesasa not installed. Run: pip install freesasa")
        return
    except Exception as e:
        print(f"[ERROR] FreeSASA analysis failed: {e}")
        return

    # 4. FoldX ddG (real, requires configured binary path)
    try:
        ddg = run_foldx(repaired_pdb_path, archive_dir)
        print(f"[SUCCESS] FoldX ddG: {ddg:.2f} kcal/mol")
    except Exception as e:
        print(f"[ERROR] FoldX failed: {e}")
        print("[INFO] Check FOLDX_BINARY at the top of this file, and that you've completed")
        print("[INFO] the free academic registration/download at foldxsuite.crg.eu")
        return

    # 5. APBS electrostatic map (real, requires configured binary path)
    try:
        apbs_map_path = run_apbs(protonated_pqr_path, archive_dir)
        print(f"[SUCCESS] APBS map: {os.path.relpath(apbs_map_path, output_base)}")
    except Exception as e:
        print(f"[ERROR] APBS failed: {e}")
        print("[INFO] Check APBS_BINARY at the top of this file, and that APBS is installed")
        print("[INFO] (conda install -c bioconda -c conda-forge apbs)")
        return

    # 6. Aggrescan3D, CamSol (structural), DeepSol -- manual (no free local/API option)
    manual_results_path = os.path.join(archive_dir, f"{safe_name}_manual_results.json")
    manual_results = load_manual_results(manual_results_path)
    if manual_results is None:
        run_manual_prepare(sequence, repaired_pdb_path, manual_results_path)
        return

    print(f"[SUCCESS] Manual results loaded from: {os.path.relpath(manual_results_path, output_base)}")
    print("-" * 110)

    results = evaluate_structural_solubility(sequence, repaired_pdb_path, sasa_result, ddg, manual_results)

    print(f"{'METRIC':<30} | {'VALUE':<15} | {'THRESHOLD / TARGET'}")
    print("-" * 110)
    print(f"{'DeepSol Probability':<30} | {results['DeepSol_Prob']:<15.2f} | >= 0.50")
    print(f"{'Aggrescan3D Patch (pLDDT>=70)':<30} | {results['Agg_Patch_Len']:<15} | <= 8 aa")
    print(f"{'CamSol (structural) Patch':<30} | {results['CamSol_Patch_Len']:<15} | <= 8 aa")
    print(f"{'Hydrophobic Fraction':<30} | {results['Hydrophobic_Fraction']:<15.2f} | < 0.25")
    print(f"{'Exposed Hydro Area':<30} | {results['Exposed_Area']:<15.1f} | < 250 \u00c5\u00b2")
    ddg_label = "FoldX \u0394\u0394G"
    print(f"{ddg_label:<30} | {results['FoldX_ddG']:<15.2f} | < +1.5 kcal/mol")
    print("-" * 110)
    print(f"FINAL DECISION : [{results['Status']}] - {results['Reason_Code']}")
    print("-" * 110)

    ts = datetime.now().strftime("%Y%m%d_%H%M")
    report_path = os.path.join(output_base, f"Step2C_Solubility_Report_{ts}.csv")

    results["Repaired_PDB_File"] = os.path.relpath(repaired_pdb_path, output_base)
    results["Protonated_PQR_File"] = os.path.relpath(protonated_pqr_path, output_base)
    results["APBS_Map_File"] = os.path.relpath(apbs_map_path, output_base)
    results["Manual_Results_File"] = os.path.relpath(manual_results_path, output_base)

    csv_data = {"Variant": winner_name, **results}
    with open(report_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=csv_data.keys())
        w.writeheader()
        w.writerow(csv_data)

    total_time = common.format_time(time.time() - start_time)
    common.print_banner("STEP 2C COMPLETE")
    print("[SUCCESS] Structural solubility & aggregation propensity checked with real tools.")
    print(f"[SUCCESS] Execution Time : {total_time}")
    print(f"[INFO] Report Saved      : {os.path.relpath(report_path, project_root)}")
    print("=" * 110 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 2C: 3D Solubility & Structural Integrity")
    parser.add_argument("--manual-prepare", action="store_true",
                         help="Print manual submission instructions and write the results template early")
    args = parser.parse_args()

    if args.manual_prepare:
        project_root = _PROJECT_ROOT
        input_csv_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
        variant_fasta_path = common.phase1g_fasta_path(project_root)
        winner_row, _ = common.get_winner_from_filtered_csv(input_csv_dir)
        if winner_row:
            variants = common.load_multi_fasta(variant_fasta_path)
            sequence = common.lookup_sequence(variants, winner_row["Variant"]) if variants else None
            safe_name = common.sanitize_variant_name(winner_row["Variant"])
            archive_dir = os.path.join(project_root, "Step_Outputs", "Phase2", "StepC", "Supplementary_Archive")
            repaired_pdb_path = os.path.join(archive_dir, f"{safe_name}_minimized.pdb")
            manual_results_path = os.path.join(archive_dir, f"{safe_name}_manual_results.json")
            if sequence:
                run_manual_prepare(sequence, repaired_pdb_path, manual_results_path)
    else:
        run_step2c_solubility_analysis()
