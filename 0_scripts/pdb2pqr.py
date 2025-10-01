import os
import sys

def run_pdb2pqr(pdb_file, pH):
    """Run PDB2PQR on a single .pdb file with a given pH and modify the resulting APBS input file."""
    if not pdb_file.endswith('.pdb'):
        return f"Skipped non-PDB file: {pdb_file}"

    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    output_pqr = f"{pdb_id}.pqr"
    output_in = f"{pdb_id}.in"
    log_file = f"{pdb_id}_pdb2pqr.log"

    # Clear existing log file
    if os.path.exists(log_file):
        os.remove(log_file)

    # Build and execute PDB2PQR command
    cmd = (
        f"stdbuf -oL pdb2pqr --ff=AMBER --whitespace --titration-state-method=propka "
        f"--with-ph={pH} --apbs-input {output_in} {pdb_file} {output_pqr} >> {log_file} 2>&1"
    )
    os.system(cmd)

    # Verify output files
    if not (os.path.exists(output_pqr) and os.path.exists(output_in)):
        raise RuntimeError(f"pdb2pqr failed for {pdb_file}. Check {log_file}")

    # Modify APBS input
    MODIFICATIONS = {
        "ion_settings": [
            "\tion charge 1 conc 0.15 radius 2.0\n",
            "\tion charge -1 conc 0.15 radius 1.8\n"
        ],
        "replacement": ("pdie 2.0000", "\tpdie 6.0\n"),
    }

    with open(output_in, 'r') as f:
        content = f.readlines()

    new_content = []
    fgcent_found = False
    for line in content:
        # Replace dielectric constant
        if line.strip() == MODIFICATIONS["replacement"][0]:
            new_content.append(MODIFICATIONS["replacement"][1])
        # Update potential output naming
        elif line.strip().startswith("write pot dx"):
            new_content.append(f"    write pot dx {pdb_id}-pot\n")
        else:
            new_content.append(line)

        # Add ion parameters after fgcent line
        if "fgcent mol 1" in line and not fgcent_found:
            new_content.extend(MODIFICATIONS["ion_settings"])
            fgcent_found = True

    with open(output_in, 'w') as f:
        f.writelines(new_content)

    print(f"Processed and customized {pdb_id} at pH {pH}")

# ============================== #
#               Main             #
# ============================== #

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python run_pdb2pqr.py <scripts_dir> <pdb_file> <pH>")
        sys.exit(1)
    
    scripts_dir = sys.argv[1]
    pdb_file = sys.argv[2]
    pH = float(sys.argv[3])

    sys.path.append(scripts_dir)
    run_pdb2pqr(pdb_file, pH)
