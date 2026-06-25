import os
import glob
import re

# directory where input .json files live
INPUT_DIR = "./"          # adjust if not current directory
RUNS_DIR = "runs"         # main runs folder
OUTPUT_LIST = "missing_runs.txt"
ENERGY_RESULTS = "ccsd_energies.txt"

# regex to match the CCSD energy line
energy_pattern = re.compile(r"CCSD total energy / hartree\s*=\s*(-?\d+\.\d+)")

missing = []
energies = []

# loop over all .json input files
for json_file in glob.glob(os.path.join(INPUT_DIR, "*.json")):
    filename = os.path.basename(json_file)
    name_no_ext = os.path.splitext(filename)[0]  # e.g. para_theta_93.91_phi_93.91
    run_folder = os.path.join(RUNS_DIR, f"run_{name_no_ext}")

    if os.path.isdir(run_folder):
        # find tamm_test.* file (should only be one, but we'll glob)
        outputs = glob.glob(os.path.join(run_folder, "tamm_test.*"))
        if not outputs:
            print(f"No tamm_test file found in {run_folder}")
            missing.append(filename)
            continue

        # take the first match
        output_file = outputs[0]
        with open(output_file, "r") as f:
            content = f.read()

        match = energy_pattern.search(content)
        if match:
            energy_val = match.group(1)
            energies.append(f"{filename}\t{energy_val}")
        else:
            print(f"No CCSD energy found in {output_file}")
    else:
        missing.append(filename)

# write missing file list
with open(OUTPUT_LIST, "w") as f:
    f.write("\n".join(missing))

# write energy results
with open(ENERGY_RESULTS, "w") as f:
    f.write("\n".join(energies))

print(f"Done. Missing files written to {OUTPUT_LIST}, energies written to {ENERGY_RESULTS}")

