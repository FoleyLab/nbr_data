import os
import sys
import subprocess

base_dir = os.path.abspath(sys.argv[1])
print(base_dir)
mdir = base_dir + "/runs"
print(mdir)

input_files = [os.path.abspath(os.path.join(base_dir, file)) for file in os.listdir(base_dir) if file.endswith('.json')]
cores = 4 

print(input_files)

for ifile in input_files:
    mname = os.path.splitext(os.path.basename(ifile))[0]
    directory_name = f"{mdir}/run_{mname}"
    print("directory_name:", directory_name)
    os.makedirs(directory_name, exist_ok=True)
    os.chdir(directory_name)
    jobscript = base_dir + "/run_script.sh"
    command = f"sbatch {jobscript} {cores} '{ifile}' {directory_name}"
    print("Executing command:", command)
    subprocess.run(command, shell=True)
