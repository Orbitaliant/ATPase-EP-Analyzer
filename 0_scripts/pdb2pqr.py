import os
import sys
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_single_pdb(pdb_file, input_directory, output_directory, log_directory):
    """Process a single PDB file with PDB2PQR."""
    if not pdb_file.endswith('.pdb'):
        return f"Skipped non-PDB file: {pdb_file}"

    pdb_id = os.path.splitext(pdb_file)[0]
    pdb_output_dir = os.path.join(output_directory, pdb_id)
    os.makedirs(pdb_output_dir, exist_ok=True)

    log_file = os.path.join(log_directory, f"{pdb_id}_pdb2pqr.log")
    input_path = os.path.join(input_directory, pdb_file)
    output_path = os.path.join(pdb_output_dir, f"{pdb_id}.pqr")
    apbs_path = os.path.join(pdb_output_dir, f"{pdb_id}.in")

    # Clear existing log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # Build and execute PDB2PQR command
    cmd = (
        f"stdbuf -oL pdb2pqr --ff=AMBER --whitespace --titration-state-method=propka "
        f"--apbs-input {apbs_path} {input_path} {output_path} >> {log_file} 2>&1"
    )
    os.system(cmd)

    # Verify output files
    if not (os.path.exists(output_path) and os.path.exists(apbs_path)):
        raise RuntimeError(f"Processing failed for {pdb_id}. Check {log_file}")

    return f"Processed {pdb_id}"

def pdb2pqr_parallel(input_dir, output_dir, log_dir):
    """Process PDB files in parallel using PDB2PQR with dynamic max_workers."""
    # Create base directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Identify PDB files
    pdb_files = [f for f in os.listdir(input_dir) if f.endswith(".pdb")]
    if not pdb_files:
        print("No PDB files found in input directory.")
        return

    # Determine max_workers dynamically
    num_pdb_files = len(pdb_files)
    num_cpu_cores = multiprocessing.cpu_count()
    max_workers = min(num_pdb_files, num_cpu_cores)

    print(
        f"\nParallel PDB2PQR processing initiated\n"
        f"- Input Directory: {input_dir}\n"
        f"- CPU Cores Used: {max_workers}\n"
        f"Note: You can adjust the number of CPU cores used by modifying `max_workers` in the script for optimum performance."
    )

    # Parallel execution
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_single_pdb, f, input_dir, output_dir, log_dir): f
            for f in pdb_files
        }

        # Track completion
        for future in as_completed(futures):
            file = futures[future]
            try:
                result = future.result()
                print(result)
            except Exception as e:
                print(f"ERROR processing {file}: {str(e)}")

    print("\nAll PDB files processed with PDB2PQR.")

# Execute parallel processing with dynamic worker allocation
pdb2pqr_parallel('3_fixed_pdb', '4_pqr_apbs_files', 'pdb2pqr_log_files')

# Modify APBS input files
apbs_output_dir = '4_pqr_apbs_files'
MODIFICATIONS = {
    "ion_settings": [
        "\tion charge 1 conc 0.15 radius 2.0\n",
        "\tion charge -1 conc 0.15 radius 1.8\n"
    ],
    "replacement": ("pdie 2.0000", "\tpdie 6.0\n"),
}

for pdb_dir in os.listdir(apbs_output_dir):
    dir_path = os.path.join(apbs_output_dir, pdb_dir)
    if not os.path.isdir(dir_path):
        continue

    for file in os.listdir(dir_path):
        if not file.endswith(".in"):
            continue

        file_path = os.path.join(dir_path, file)
        with open(file_path, 'r') as f:
            content = f.readlines()

        new_content = []
        fgcent_found = False
        for line in content:
            # Replace dielectric constant
            if line.strip() == MODIFICATIONS["replacement"][0]:
                new_content.append(MODIFICATIONS["replacement"][1])
            # Update potential output naming
            elif line.strip().startswith("write pot dx"):
                new_content.append(f"    write pot dx {pdb_dir}-pot\n")
            else:
                new_content.append(line)

            # Add ion parameters after fgcent line
            if "fgcent mol 1" in line and not fgcent_found:
                new_content.extend(MODIFICATIONS["ion_settings"])
                fgcent_found = True

        with open(file_path, 'w') as f:
            f.writelines(new_content)

print("\nAPBS input files customized successfully.")
