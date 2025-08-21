import os
import multiprocessing
import MDAnalysis as mda
import pandas as pd
import gc

def process_pdb_entry(entry):
    """Processes a single PDB entry, extracting chains and writing output."""
    accession, chain_list, target_dir = entry
    file_path = get_pdb(accession)
    u = None
    if not os.path.exists(file_path):
        print(f"Skipping {file_path}: File not found")
        return False

    try:
        # load only once
        u = mda.Universe(file_path)
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return False
    try:
        for chain_ID in chain_list:
            output_pdb = os.path.join(target_dir, f"{accession}_{chain_ID}.pdb")
            if os.path.exists(output_pdb):
                print(f"Skipping {output_pdb}: already exists")
                continue
            try:
                chain = u.select_atoms(f"protein and chainID {chain_ID}")
                if chain.n_atoms > 0:
                    #output_pdb = os.path.join(target_dir, f"{accession}_{chain_ID}.pdb")
                    with mda.Writer(output_pdb) as w:
                        w.write(chain)
                    post_process(output_pdb)
                    print(f"Processed {output_pdb}")
                else:
                    print(f"Chain {chain_ID} not found in {file_path}")
            except Exception as e:
                print(f"Error processing {file_path} (Chain {chain_ID}): {e}")
    # release memory
    finally:
        del u
        gc.collect()
    
    return True

def strip_to_chain(pdb_file, chain_ID):
    """Depracated, has been re-written into `process_pdb_entry`"""
    try:
        u = mda.Universe(pdb_file)
        chain = u.select_atoms(f"protein and chainID {chain_ID}")
        if len(chain) == 0:
            return None
        return chain
    except Exception as e:
        print(f"Error loading {pdb_file}: {e}")
        return None

def post_process(fname):
    """Removes unnecessary 'TER' lines from the PDB file to clean the output"""
    try:
        with open(fname, "r") as f:
            lines = f.readlines()

        # Keep only necessary TER lines (last two lines)
        final_lines = lines[-2:]
        no_ter = [line for line in lines if not line.startswith("TER") or line in final_lines]

        if len(no_ter) != len(lines):
            with open(fname, "w") as f:
                f.writelines(no_ter)
    except Exception as e:
        print(f"Error post-processing {fname}: {e}")

def get_pdb(id, path='./PDBs/'):
    """Generate a cross-platform PDB file path"""
    return os.path.join(path, f"{id}.pdb")

def parallel_pdb_processing(pdb_data, target_dir):
    """Processes multiple PDB files in parallel"""
    os.makedirs(target_dir, exist_ok=True)

    # Filter downloaded PDB entries
    pdb_entries = pdb_data[pdb_data['Downloaded'] == True].copy()

    # Convert chain list from semicolon-separated string to a list
    pdb_entries['Chain_list'] = pdb_entries['Chains'].apply(lambda x: x.split(';'))

    # Prepare input list for multiprocessing
    task_list = [(row['Accession'], row['Chain_list'], target_dir) for _, row in pdb_entries.iterrows()]

    num_workers = min(multiprocessing.cpu_count(), 8)

    with multiprocessing.Pool(processes=num_workers, maxtasksperchild=10) as pool:
        pool.map(process_pdb_entry, task_list)

    print("All PDB processing completed.")

if __name__ == "__main__":
    target_dir = "Results/activation_segments/unaligned"
    tsv_path = 'structure-matching-IPR011009.tsv'
    folder_path = './PDBs'

    # Ensure pdb_data is properly initialized before calling the function
    try:
        # Load PDB data from TSV file
        pdb_data = pd.read_csv(tsv_path, sep="\t", header=0, engine='python')
        pdb_data['Accession'] = pdb_data['Accession'].str.upper()

        # Get downloaded PDBs
        file_names = [os.path.splitext(f)[0] for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
        
        # Convert to uppercase for correct matching
        pdb_raw = pd.DataFrame({"PDBs": file_names})
        pdb_raw['PDBs'] = pdb_raw['PDBs'].str.upper()

        # Mark downloaded PDBs
        pdb_data['Downloaded'] = pdb_data['Accession'].isin(pdb_raw['PDBs'])
        parallel_pdb_processing(pdb_data, target_dir)
    except Exception as e:
        print(f"Error loading pdb_data: {e}")