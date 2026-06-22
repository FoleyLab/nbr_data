import os
import shutil
import numpy as np
import psi4
from concurrent.futures import ProcessPoolExecutor
from cqed_scf.calculator import CQEDCalculator

# --- Configuration ---
NUM_CORES = 16
MEM_PER_CORE = "4 GB"
OUTPUT_FILE = "ortho_meta_refinement_scan.txt"

# ----------------------------
# Molecular Geometries
# ----------------------------
GEOMS = {
    "ortho": """
1 1
C           -1.804928163307     1.957993763262     0.703312273806
C           -0.379708783307     1.994122833262     0.698532703806
C            0.296125016693     0.817793533262     0.710271493806
C           -2.520286433307     0.755089873262     0.736288843806
H           -2.344947113307     2.899196893262     0.691895063806
H            0.158564066693     2.933869823262     0.699142733806
H           -3.601954283307     0.764862203262     0.746931053806
N            1.767881836693     0.820900013262     0.771891313806
O            2.315054046693    -0.296733496738     0.879853723806
O            2.340645916693     1.923356243262     0.711986073806
C           -1.829967733307    -0.442167236738     0.756258983806
H           -2.356763623307    -1.389967436738     0.789740873806
C           -0.361341153307    -0.491572936738     0.714148383806
H            0.119338216693    -1.238105076738     1.350400383806
BR          -0.151212663307    -1.224162306738    -1.170925976194
units angstrom
no_reorient
no_com
symmetry c1
""",
    "meta": """
1 1
         C           -0.929257263947     2.021527608578     0.744707683350
         C            0.476075706053     1.968481358578     0.682883583350
         C            1.153033166053     0.732862858578     0.67108907335
         C            0.486309286053    -0.455398891422     0.707696283350
         C           -1.646688783947     0.850023888578     0.786483593350
         H           -1.430027043947     2.980198348578     0.754644003350
         H            1.068570756053     2.878318968578     0.644324213350
         H            1.030908186053    -1.394630481422     0.699715393350
         H           -2.730391873947     0.862207158578     0.834726773350
         N            2.627601876053     0.732774608578     0.609077593350
         O            3.188360516053    -0.377859281422     0.588451963350
         O            3.186221516053     1.845711198578     0.586422223350
         C           -0.982368843947    -0.464026221422     0.760065283350
         H           -1.395507033947    -1.190671951422     1.465426213350
         BR          -1.494673453947    -1.187920261422    -1.064256256650
units angstrom
no_reorient
no_com
symmetry c1
"""
}

psi4_options = {
    "basis": "6-311G*",
    "scf_type": "df",
    "e_convergence": 1e-10,
    "d_convergence": 1e-10,
    "dft_radial_points": 99,
    "dft_spherical_points": 590,
}

def run_point(task):
    theta, phi, worker_id = task
    scratch_path = f"/tmp/psi_work_om_{worker_id}"
    os.makedirs(scratch_path, exist_ok=True)

    psi4.core.IOManager.shared_object().set_default_path(scratch_path)
    psi4.set_memory(MEM_PER_CORE)
    psi4.core.set_num_threads(1)
    psi4.core.set_output_file(os.path.join(scratch_path, "output.dat"), False)

    theta_rad, phi_rad = np.radians(theta), np.radians(phi)
    field_vec = np.array([
        np.sin(theta_rad) * np.cos(phi_rad),
        np.sin(theta_rad) * np.sin(phi_rad),
        np.cos(theta_rad)
    ]) * 0.1

    results = {"theta": theta, "phi": phi, "vec": field_vec, "energies": {}}

    for name, geom_string in GEOMS.items():
        try:
            calc = CQEDCalculator(
                lambda_vector=field_vec,
                psi4_options=psi4_options,
                omega=0.06615,
                charge=1,
                multiplicity=1,
                functional="wb97x",
            )
            results["energies"][name] = calc.energy(geom_string)
            psi4.core.clean()
        except Exception as e:
            results["energies"][name] = np.nan
            print(f"Error at {theta}, {phi} for {name}: {e}")

    shutil.rmtree(scratch_path, ignore_errors=True)
    return results

if __name__ == "__main__":
    # Target: 73, 36 over +/- 10 degrees with 1-degree intervals (21 points per dimension)
    theta_list = np.linspace(63, 83, 21)
    phi_list = np.linspace(26, 46, 21)

    tasks = []
    counter = 0
    for t in theta_list:
        for p in phi_list:
            tasks.append((t, p, counter))
            counter += 1

    print(f"Submitting {len(tasks)} points to {NUM_CORES} cores...")

    with ProcessPoolExecutor(max_workers=NUM_CORES) as executor:
        all_data = list(executor.map(run_point, tasks))

    all_data.sort(key=lambda x: (x["theta"], x["phi"]))

    with open(OUTPUT_FILE, "w") as f:
        header = f"{'Theta':>10} {'Phi':>10} {'Ex':>12} {'Ey':>12} {'Ez':>12} {'Ortho_E':>20} {'Meta_E':>20}\n"
        f.write(header)
        f.write("-" * len(header) + "\n")

        for d in all_data:
            t, p = d["theta"], d["phi"]
            vx, vy, vz = d["vec"]
            e_o = d["energies"].get("ortho", np.nan)
            e_m = d["energies"].get("meta", np.nan)

            line = f"{t:10.3f} {p:10.3f} {vx:12.6f} {vy:12.6f} {vz:12.6f} {e_o:20.12f} {e_m:20.12f}\n"
            f.write(line)

    print(f"Scan complete. Data written to {OUTPUT_FILE}")
