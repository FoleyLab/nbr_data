import csv
import numpy as np
import psi4
from cqed_scf import CQEDCalculator
from cqed_scf.utils import write_xyz, ANGSTROM_TO_BOHR, generate_field_vector_from_theta_and_phi
# =========================
# Psi4 geometry (angstrom)
# =========================
para_string = """
1 1
C     -0.5116182968     1.2443860245     0.7321400487
C      0.8565005932     1.2519037145     0.7179482187
C      1.5241187232     0.0246619245     0.7139277887
H     -1.0718043968     2.1726823145     0.7459257087
H      1.4361289632     2.1639218745     0.7120990087
N      3.0085395832     0.0460971045     0.6988237987
O      3.5750973032    -1.0827681655     0.6991747087
O      3.5421143632     1.1908708545     0.6892020187
C     -0.4754649468    -1.2534027655     0.7421186387
H     -1.0085744268    -2.1973779555     0.7629457887
C      0.8922277032    -1.2214078055     0.7280658187
H      1.4980484232    -2.1162446955     0.7296532087
C     -1.2679065768    -0.0158418055     0.7121274187
H     -2.1165207968    -0.0252932155     1.4031614987
Br    -2.1149669868    -0.0346669255    -1.1218740813
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
theta = 101  # degrees from z-axis
phi = 31  # degrees from x-axis in xy-plane
field_vector_center = generate_field_vector_from_theta_and_phi(theta, phi)
print(F"Field Magnitude before scaling is {np.linalg.norm(field_vector_center)}")
print(F"Field Vector is")
print(field_vector_center)
lambda_direction = np.asarray(field_vector_center, dtype=float)
lambda_direction /= np.linalg.norm(lambda_direction)
omega = 0.06615

# lambda magnitudes to scan over, same (theta, phi) direction for all points
lambda_magnitude = 0.1 

HARTREE_TO_KCAL_MOL = 627.5094740631

lambda_vector = (lam_mag * lambda_direction).tolist()
Ex, Ey, Ez = lambda_vector
print(F"--- lambda magnitude = {lam_mag} ---")
print(F"Mag after scaling {np.linalg.norm(lambda_vector)}")

calc = CQEDCalculator(lambda_vector=lambda_vector,psi4_options=psi4_options,omega=omega,density_fitting=True,charge=1,multiplicity=1,functional="wb97x")

energy_para = calc.energy(para_string)
