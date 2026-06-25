import os
import glob
import re
import shutil

# === USER SETTINGS ===
HOMEDIR = os.path.expanduser("/global/cfs/projectdirs/m4262/mapol/exachem-dev/inputs/ci/NB_QED_RHF_INPUTS/SUB_META_PARA/ortho_cc_files")
INPUT_DIR = ""                       # where your .json input files live
RUN_DIRS = [                           # list all run directories to check
    os.path.join(HOMEDIR, "runs"),
    os.path.join(HOMEDIR, "needs_running", "runs")
]
#RUN_DIRS = os.path.join(HOMEDIR, "runs")
NEEDS_RUNNING_DIR = os.path.join(HOMEDIR, "needs_running")
OUTPUT_LIST = "missing_runs.txt"
ENERGY_RESULTS = "ccsd_energies.tsv"

# === PREP ===
os.makedirs(NEEDS_RUNNING_DIR, exist_ok=True)

# regex patterns
energy_pattern = re.compile(r"CCSD total energy / hartree\s*=\s*(-?\d+\.\d+(?:[Ee][+-]?\d+)?)")
theta_phi_pattern = re.compile(r"theta_([-\d\.]+)_phi_([-\d\.]+)")

missing = []
rows = []

# === MAIN LOOP ===
for json_file in glob.glob(os.path.join(INPUT_DIR, "*.json")):
    filename = os.path.basename(json_file)
    name_no_ext = os.path.splitext(filename)[0]

    # extract theta and phi
    tp_match = theta_phi_pattern.search(name_no_ext)
    if tp_match:
        theta_val, phi_val = tp_match.groups()
    else:
        theta_val, phi_val = ("NA", "NA")

    found_output = False

    # search through all possible run directories
    for base_run_dir in RUN_DIRS:
        run_folder = os.path.join(base_run_dir, f"run_{name_no_ext}")
        print(F"GOING INTO RUN FOLDER {run_folder}")
        if os.path.isdir(run_folder):
            outputs = glob.glob(os.path.join(run_folder, "tamm_test.*"))
            print(F"FOUND OUTPUT {outputs}")
            if not outputs:
                continue  # try next base dir

            # choose newest output if multiple
            output_file = max(outputs, key=os.path.getmtime)
            with open(output_file, "r") as f:
                content = f.read()

            match = energy_pattern.search(content)
            if match:
                energy_val = match.group(1)
                rows.append(f"{theta_val}\t{phi_val}\t{energy_val}")
            else:
                rows.append(f"{theta_val}\t{phi_val}\tMISSING_ENERGY")

            found_output = True
            break  # stop checking other run dirs once found

    # if not found in any directory
    if not found_output:
        missing.append(filename)
        dest = os.path.join(NEEDS_RUNNING_DIR, filename)
        if not os.path.exists(dest):  # avoid overwriting existing file
            shutil.move(json_file, dest)

# === WRITE OUTPUTS ===
with open(OUTPUT_LIST, "w") as f:
    f.write("\n".join(missing))

with open(ENERGY_RESULTS, "w") as f:
    f.write("theta\tphi\tenergy\n")
    f.write("\n".join(rows))

print(f"Done.")
print(f"- Missing files written to {OUTPUT_LIST}")
print(f"- Results written to {ENERGY_RESULTS}")
print(f"- Missing .json files moved to {NEEDS_RUNNING_DIR}/")
