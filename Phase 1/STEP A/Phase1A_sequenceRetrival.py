import os
import sys
import time
from datetime import datetime
from Bio import Entrez
from Bio import SeqIO

# =============================================================================
# EXPERIMENTAL CONFIGURATION
# =============================================================================

# NCBI Identification (Required for API compliance)
# Uses an environment variable for safety/portability, falling back to your email
Entrez.email = os.getenv("NCBI_EMAIL", "enzoleonor.3309@gmail.com") 

# Target viral proteomes for Chimeric Vaccine construct
# Queries are strictly refined for Monkeypox Clade IIa/IIb and HIV-1 CRF01_AE
TARGET_CLUSTERS = {
    "Mpox_L1R": 'L1R AND "Monkeypox virus"[Organism] AND ("clade IIa" OR "clade IIb")',
    "Mpox_B5R": 'B5R AND "Monkeypox virus"[Organism] AND ("clade IIa" OR "clade IIb")',
    "Mpox_A35R": 'A35R AND "Monkeypox virus"[Organism] AND ("clade IIa" OR "clade IIb")',
    "HIV_gp120": 'gp120 AND "HIV-1"[Organism] AND CRF01_AE',
    "HIV_gp41": 'gp41 AND "HIV-1"[Organism] AND CRF01_AE',
    "HIV_p24": 'p24 AND "HIV-1"[Organism] AND CRF01_AE',
    "HIV_p17": 'p17 AND "HIV-1"[Organism] AND CRF01_AE'
}

# Sampling parameters
VARIANTS_PER_TARGET = 30 
TOTAL_PROJECTED = len(TARGET_CLUSTERS) * VARIANTS_PER_TARGET

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def format_time(seconds):
    """Formats time in seconds to a readable MM:SS format."""
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def run_high_density_retrieval():
    start_time = time.time()
    
    # 1. DYNAMIC DIRECTORY RESOLUTION
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    phase1a_path = os.path.join(project_root, "Step_Outputs", "Phase1A")
    
    os.makedirs(phase1a_path, exist_ok=True)

    # 2. EXPERIMENTAL PREVIEW & HEADER
    print("\n" + "="*80)
    print(f"{'PHASE 1A: VIRAL PROTEOME SEQUENCE RETRIEVAL (NCBI ENTREZ)':^80}")
    print("="*80)
    print(f"[INFO] Initialization Time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Target Antigens     : {len(TARGET_CLUSTERS)} defined")
    print(f"[INFO] Variants per Target : {VARIANTS_PER_TARGET}")
    print(f"[INFO] Projected Yield     : {TOTAL_PROJECTED} FASTA sequences")
    print(f"[INFO] Output Directory    : {phase1a_path}")
    print("-" * 80)

    # 3. WORKSPACE PREPARATION (Purging old data)
    sys.stdout.write("[PROCESS] Purging prior sequence data to ensure experimental integrity...")
    sys.stdout.flush()
    purged_count = 0
    for f in os.listdir(phase1a_path):
        if f.endswith(".fasta"):
            os.remove(os.path.join(phase1a_path, f))
            purged_count += 1
    print(f" Done. ({purged_count} files removed)")
    print("-" * 80)

    # 4. SEQUENCE ACQUISITION PROTOCOL
    successful_downloads = 0

    for label, query in TARGET_CLUSTERS.items():
        print(f"\n[INFO] Establishing NCBI connection for target: {label}")
        
        try:
            # 4a. Query execution
            search_handle = Entrez.esearch(db="protein", term=query, retmax=VARIANTS_PER_TARGET)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            id_list = search_results["IdList"]
            retrieved_count = len(id_list)
            
            if retrieved_count == 0:
                print(f"[WARNING] No records found for query: '{query}'")
                continue

            # 4b. Batch Data Fetching
            elapsed_search = format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ FETCH ] {label} | Found {retrieved_count:02d} IDs | Downloading batch... | Elapsed: {elapsed_search}")
            sys.stdout.flush()

            # Execute a single batch request for all IDs
            fetch_handle = Entrez.efetch(db="protein", id=",".join(id_list), rettype="fasta", retmode="text")
            records = list(SeqIO.parse(fetch_handle, "fasta"))
            fetch_handle.close()

            # 4c. Parsing and saving files individually
            for i, record in enumerate(records):
                file_name = f"{label}_Var_{i+1:02d}_{record.id}.fasta"
                with open(os.path.join(phase1a_path, file_name), "w") as f:
                    SeqIO.write(record, f, "fasta")
                successful_downloads += 1

            elapsed_done = format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ FETCH ] {label} | Successfully saved {len(records):02d} records | Elapsed: {elapsed_done:<15}\n")
            sys.stdout.flush()

            # Regulated delay to comply with NCBI API restrictions between cluster requests
            time.sleep(1.0) 

        except Exception as e:
            print(f"\n[ERROR] Protocol failure during {label} acquisition. Reason: {e}")

    # 5. POST-EXECUTION REPORT
    total_time = format_time(time.time() - start_time)
    print("\n" + "="*80)
    print(f"{'ACQUISITION PROTOCOL COMPLETE':^80}")
    print("="*80)
    print(f"[SUCCESS] Total Sequences Retrieved : {successful_downloads}/{TOTAL_PROJECTED}")
    print(f"[SUCCESS] Total Execution Time      : {total_time}")
    print(f"[INFO] Data formatting complete. Proceed to Phase 1B for stability thresholding.")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_high_density_retrieval()
