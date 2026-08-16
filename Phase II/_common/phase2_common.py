"""
phase2_common.py

Shared utilities for the Phase II vaccine-candidate screening pipeline
(Steps 2A-2D).

WHY THIS FILE EXISTS:
Every bug found while debugging Steps 2A-2D this session came from the
same root cause: path-resolution and FASTA-parsing logic was copy-pasted
into each script, and the copies drifted out of sync with each other
(one script got fixed, the others didn't). Centralizing that logic here
means it only needs to be correct once. If you need to change how the
project root is found, or how FASTA headers are parsed, change it here
and all four steps pick it up automatically.

Each Step_2X script keeps a small (~10 line) local bootstrap that finds
this file and adds it to sys.path -- that bootstrap is intentionally
kept in each script rather than moved here, since a module can't be
used to locate itself.
"""

import os
import sys
import csv

RESEARCH_ANCHOR = "Research"
PHASE1G_FASTA_RELATIVE = os.path.join(
    "Phase 1", "Step_Outputs", "Phase1G", "Phase1G_Constructs_2026-07-22_2251.fasta"
)


def resolve_project_root(script_file):
    """
    Walk upward from a script's own location until a folder literally
    named "Research" is found. Anchoring on a named folder (rather than
    a fixed hop count) means this keeps working even though STEP A/B/C/D
    scripts don't all sit at the same folder depth -- which is exactly
    what broke the fixed "../../.." version of this logic.
    """
    script_dir = os.path.dirname(os.path.abspath(script_file))
    current = script_dir
    while os.path.basename(current) != RESEARCH_ANCHOR:
        parent = os.path.dirname(current)
        if parent == current:
            print(f"\n[FATAL ERROR] Could not locate a '{RESEARCH_ANCHOR}' anchor folder above: {script_dir}")
            sys.exit(1)
        current = parent
    return current


def phase1g_fasta_path(project_root):
    """Absolute path to the single multi-sequence Phase 1G output FASTA."""
    return os.path.join(project_root, PHASE1G_FASTA_RELATIVE)


def load_multi_fasta(fasta_path):
    """
    Parses a multi-sequence FASTA file into an ordered dict of
    {header_text: sequence}. Header text is everything after '>' on the
    header line, stripped of whitespace and otherwise left untouched
    (metadata like "| length=138" is preserved as-is so it matches
    exactly whatever Step 2A wrote into its CSV's Variant column).

    Returns None if the file doesn't exist (caller decides how to
    report that), or an (possibly empty) dict otherwise.
    """
    if not os.path.isfile(fasta_path):
        return None

    records = {}
    current_name = None
    current_seq_lines = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records[current_name] = "".join(current_seq_lines).upper()
                current_name = line[1:].strip()
                current_seq_lines = []
            else:
                current_seq_lines.append(line)
        if current_name is not None:
            records[current_name] = "".join(current_seq_lines).upper()

    return records


def sanitize_variant_name(name):
    """
    Variant names carry '| length=NNN' metadata and spaces (e.g.
    "Vax_Var1_fb3faeb6 | length=138"), which is fragile to put directly
    into filenames. This strips the metadata and swaps spaces for
    underscores, producing one stable filename-safe token that Steps
    2B, 2C, and 2D all build filenames from identically.
    """
    return name.split("|")[0].strip().replace(" ", "_")


def lookup_sequence(variants, full_name):
    """
    Looks up a sequence by trying the exact recorded name first, then
    falling back to the name with any '| metadata' stripped, in case
    the FASTA header and a CSV's Variant column ever drift slightly.
    """
    if variants is None:
        return None
    clean_name = full_name.split("|")[0].strip()
    return variants.get(full_name) or variants.get(clean_name)


def get_winner_from_filtered_csv(filtered_dir):
    """
    Finds the most recent Filtered CSV written by Step 2A and returns
    (winner_row_dict, csv_path) for the rank-1 (first) row -- Step 2A
    always writes this file sorted best-first by Stability Index.

    Returns (None, None) after printing a diagnostic message if
    anything is missing, so every downstream step reports failures
    the same way instead of each reinventing this check slightly
    differently (which is how Steps 2B/2C/2D ended up with three
    subtly different error-handling styles before this rewrite).
    """
    if not os.path.isdir(filtered_dir):
        print(f"[ERROR] Filtered output folder does not exist: {filtered_dir}")
        print("[ERROR] Run Step 2A first so this folder gets created.")
        return None, None

    csv_files = sorted([f for f in os.listdir(filtered_dir) if f.endswith(".csv")])
    if not csv_files:
        print(f"[ERROR] No filtered candidate CSVs found in: {filtered_dir}")
        print("[ERROR] Step 2A ran but produced zero viable candidates -- check Rejection_Reasons in its Raw log.")
        return None, None

    csv_path = os.path.join(filtered_dir, csv_files[-1])
    with open(csv_path, 'r') as f:
        reader = list(csv.DictReader(f))

    if not reader:
        print(f"[ERROR] Filtered CSV is empty: {csv_path}")
        return None, None

    return reader[0], csv_path


def print_banner(text, width=115):
    print("\n" + "=" * width)
    print(f"{text:^{width}}")
    print("=" * width)


def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"


# ---------------------------------------------------------------------------
# Contextual-only physicochemical helpers not provided directly by Biopython
# ---------------------------------------------------------------------------

def compute_aliphatic_index(seq):
    """
    Ikai (1980) aliphatic index -- relative volume occupied by aliphatic
    side chains (Ala, Val, Ile, Leu). Biopython's ProteinAnalysis does
    NOT implement this, so it's computed by hand here.

    Per methodology this is a CONTEXTUAL metric only -- it must never be
    used as a rejection criterion.
    """
    length = len(seq)
    if length == 0:
        return 0.0
    ala = seq.count('A') / length * 100
    val = seq.count('V') / length * 100
    ile = seq.count('I') / length * 100
    leu = seq.count('L') / length * 100
    return ala + 2.9 * val + 3.9 * (ile + leu)


# ExPASy ProtParam's mammalian (in vitro) N-end rule half-life table,
# in hours, keyed by N-terminal residue (Bachmair et al. 1986;
# Rogers et al. 1986).
_MAMMALIAN_HALF_LIFE_HOURS = {
    'A': 4.4, 'R': 1.0, 'N': 1.4, 'D': 1.1, 'C': 1.2,
    'Q': 0.8, 'E': 1.0, 'G': 30.0, 'H': 3.5, 'I': 20.0,
    'L': 5.5, 'K': 1.3, 'M': 30.0, 'F': 1.1, 'P': 20.0,
    'S': 1.9, 'T': 7.2, 'W': 2.8, 'Y': 2.8, 'V': 100.0,
}


def estimate_half_life_hours(seq):
    """
    Rough N-end-rule half-life estimate (mammalian, in vitro) based on
    the sequence's N-terminal residue, matching the convention used by
    ExPASy ProtParam. This is a approximation, not a real assay result.

    Per methodology this is a CONTEXTUAL metric only -- it must never be
    used as a rejection criterion.
    """
    if not seq:
        return None
    return _MAMMALIAN_HALF_LIFE_HOURS.get(seq[0])
