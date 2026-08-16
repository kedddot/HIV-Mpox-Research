import os
import csv
import random
import hashlib
from datetime import datetime

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def print_banner(text): 
    print(f"\n{'='*80}\n{text:^80}\n{'='*80}")

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def run_step1g_construction():
    # 1. PATH RESOLUTION
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1F", "Filtered")
    output_dir = os.path.join(project_root, "Step_Outputs", "Phase1G")
    os.makedirs(output_dir, exist_ok=True)

    print_banner("PHASE 1G: CHIMERIC VACCINE ASSEMBLY & LINKER INTEGRATION")

    # 2. LOAD ELITE EPITOPES
    try:
        files = sorted([f for f in os.listdir(input_folder) if f.endswith(".csv")])
        latest_csv = os.path.join(input_folder, files[-1])
    except IndexError:
        print("[ERROR] No input CSV found in Phase 1F Filtered directory.")
        return

    mhc_i, mhc_ii, bcell = [], [], []
    
    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pep = row['Peptide']
            pep_type = row['Type']
            if pep_type == "MHC-I" and pep not in mhc_i:
                mhc_i.append(pep)
            elif pep_type == "MHC-II" and pep not in mhc_ii:
                mhc_ii.append(pep)
            elif pep_type == "B-cell" and pep not in bcell:
                bcell.append(pep)

    print(f"[INFO] Source Data       : {os.path.basename(latest_csv)}")
    print(f"[INFO] Available Pool    : {len(mhc_i)} MHC-I | {len(mhc_ii)} MHC-II | {len(bcell)} B-cell")
    print("-" * 80)

    # 3. VACCINE CONFIGURATION
    ADJUVANT = "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK" # Beta-defensin 3
    ADJ_LINKER = "EAAAK"
    
    # Class-specific linkers
    L_MHCI = "AAY"
    L_MHCII = "GPGPG"
    L_BCELL = "KK"

    constructs = []
    
    # Generate 5 variant constructs
    for v in range(1, 6):
        # Sample proportionally to ensure a balanced multi-epitope construct
        # (Defaults to taking max available if less than 5 are present in a class)
        sample_i = random.sample(mhc_i, min(5, len(mhc_i)))
        sample_ii = random.sample(mhc_ii, min(5, len(mhc_ii)))
        sample_b = random.sample(bcell, min(5, len(bcell)))

        # Link intra-class epitopes
        block_i = L_MHCI.join(sample_i)
        block_ii = L_MHCII.join(sample_ii)
        block_b = L_BCELL.join(sample_b)

        # Assemble the final chain. 
        # Transitions between blocks utilize the linker of the incoming block.
        chain_blocks = []
        if block_i: chain_blocks.append(block_i)
        if block_ii: chain_blocks.append(block_ii)
        if block_b: chain_blocks.append(block_b)

        # Join the major blocks safely based on what's available
        epitope_chain = ""
        if sample_i:
            epitope_chain += block_i
        if sample_ii:
            epitope_chain += (L_MHCII if epitope_chain else "") + block_ii
        if sample_b:
            epitope_chain += (L_BCELL if epitope_chain else "") + block_b

        # Final Assembly: Adjuvant + EAAAK + Epitope Chain
        final_sequence = f"{ADJUVANT}{ADJ_LINKER}{epitope_chain}"
        seq_hash = hashlib.md5(final_sequence.encode()).hexdigest()[:8]
        
        constructs.append({
            "Construct_ID": f"Vax_Var{v}_{seq_hash}",
            "Length": len(final_sequence),
            "MHC_I_Count": len(sample_i),
            "MHC_II_Count": len(sample_ii),
            "BCell_Count": len(sample_b),
            "Sequence": final_sequence
        })

        print(f"[SUCCESS] Generated Variant {v} ({seq_hash}) | Length: {len(final_sequence)} aa")

    # 4. EXPORT
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    csv_path = os.path.join(output_dir, f"Phase1G_Constructs_{ts}.csv")
    fasta_path = os.path.join(output_dir, f"Phase1G_Constructs_{ts}.fasta")

    # Save to CSV for data matrix
    with open(csv_path, 'w', newline='') as f:
        fieldnames = ["Construct_ID", "Length", "MHC_I_Count", "MHC_II_Count", "BCell_Count", "Sequence"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(constructs)

    # Save to FASTA for downstream modeling (Phase 2/3)
    with open(fasta_path, 'w') as f:
        for c in constructs:
            # Wrap sequence to 80 characters per line for standard FASTA format
            seq = c['Sequence']
            wrapped_seq = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
            f.write(f">{c['Construct_ID']} | length={c['Length']}\n{wrapped_seq}\n")

    print("-" * 80)
    print(f"[INFO] Exported Matrix : {os.path.basename(csv_path)}")
    print(f"[INFO] Exported FASTA  : {os.path.basename(fasta_path)}")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_step1g_construction()