import os
import glob
import re

# directory where input .json files live
INPUT_DIR = "./"          # adjust if not current directory
RUNS_DIR = "runs"         # main runs folder
OUTPUT_LIST = "missing_runs.txt"
ENERGY_RESULTS = "ccsd_energies.tsv"   # use tab-separated values for easy loading

# regex to match the CCSD energy line
energy_pattern = re.compile(r"CCSD total energy / hartree\s*=\s*(-?\d+\.\d+)")
# regex to extract theta and phi
theta_phi_pattern = re.compile(r"theta_([-\d\.]+)_phi_([-\d\.]+)")

missing = []
rows = []

# loop over all .json input files
for json_file in glob.glob(os.path.join(INPUT_DIR, "*.json")):
    filename = os.path.basename(json_file)
    name_no_ext = os.path.splitext(filename)[0]  # e.g. para_theta_93.91_phi_93.91
    run_folder = os.path.join(RUNS_DIR, f"run_{name_no_ext}")

    # extract theta and phi
    tp_match = theta_phi_pattern.search(name_no_ext)
    if tp_match:
        theta_val, phi_val = tp_match.groups()
    else:
        theta_val, phi_val = ("NA", "NA")

    if os.path.isdir(run_folder):
        # find tamm_test.* file
        outputs = glob.glob(os.path.join(run_folder, "tamm_test.*"))
        if not outputs:
            missing.append(filename)
            continue

        output_file = outputs[0]
        with open(output_file, "r") as f:
            content = f.read()

        match = energy_pattern.search(content)
        if match:
            energy_val = match.group(1)
            rows.append(f"{theta_val}\t{phi_val}\t{energy_val}")
        else:
            rows.append(f"{theta_val}\t{phi_val}\tMISSING_ENERGY")
    else:
        missing.append(filename)

# write missing file list
with open(OUTPUT_LIST, "w") as f:
    f.write("\n".join(missing))

# write theta/phi/energy results
with open(ENERGY_RESULTS, "w") as f:
    f.write("theta\tphi\tenergy\n")   # header row
    f.write("\n".join(rows))

print(f"Done. Missing files written to {OUTPUT_LIST}, results written to {ENERGY_RESULTS}")

