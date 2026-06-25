"""
Geometry optimization of water in an optical cavity using CQED-RHF + BFGS,
with Cartesian projected-gradient removal of translation/rotation.
"""

import numpy as np
import psi4
#psi4.core.be_quiet()

from cqed_scf import CQEDCalculator
from cqed_scf.drivers import bfgs_optimize
from cqed_scf.utils import write_xyz, ANGSTROM_TO_BOHR, generate_field_vector_from_theta_and_phi
# =========================
# Psi4 geometry (angstrom)
# =========================

meta_string = """
1 1
C           -0.929257263947     2.021527608578     0.744707683350
C            0.476075706053     1.968481358578     0.682883583350 
C            1.153033166053     0.732862858578     0.671089073350 
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
Br          -1.494673453947    -1.187920261422    -1.064256256650 
units angstrom
no_reorient
no_com
symmetry c1
"""

# =========================
# Psi4 options
# =========================

psi4_options = {
    "basis": "6-311G*",
    "scf_type": "df",
    "e_convergence": 1e-12,
    "d_convergence": 1e-12,
    "dft_radial_points": 99,
    "dft_spherical_points": 590,
    "dft_pruning_scheme": "none"
}

psi4.set_options(psi4_options)

# =========================
# Cavity parameters
# =========================

theta_central = 70 # 70° from z-axis
phi_central = 31 # 31° from x-axis in xy-plane
d_alpha = 1.0 # deviation angle in degrees

# pre-compute different field vectors for finite differences of QED-RHF energy wrt theta and phi
field_vector_center = generate_field_vector_from_theta_and_phi(theta_central, phi_central)
vec_mag = np.linalg.norm(field_vector_center)
print(F"Field Magnitude before scaling is {vec_mag}")
print(F"Field Vector is")
print(field_vector_center)

lambda_direction = np.asarray(field_vector_center, dtype=float)
lambda_direction /= np.linalg.norm(lambda_direction)
lam_mag = 0.1
omega = 0.06615

lambda_vector = lam_mag * lambda_direction

calc = CQEDCalculator(
        lambda_vector=lambda_vector,
        psi4_options=psi4_options,
        omega=omega,
        density_fitting=True,
        charge=1,
        multiplicity=1,
        functional="wb97x",  # try None for RHF
)

E = calc.energy(meta_string)

print(F"QED-DFT Energy is {E:16.12f}")
