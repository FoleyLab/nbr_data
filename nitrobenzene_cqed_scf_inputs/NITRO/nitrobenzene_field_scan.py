import numpy as np
import psi4

from cqed_rhf.utils import write_xyz
from cqed_rhf.calculator import CQEDRHFCalculator
from cqed_rhf.drivers import velocity_verlet_md
from cqed_rhf.observables.nitrobenzene_orientation import NitrobenzeneOrientation
from cqed_rhf.utils import write_xyz, ANGSTROM_TO_BOHR

#<-- Change this to a number of theta and phi values
#<-- you want to compute the scan over!
num_phi_vals = 2
num_theta_vals = 2

def generate_field_vector_from_theta_and_phi(theta, phi):
    """
    Generate a unit field vector from spherical coordinates.

    Parameters:
    -----------
    theta : float
        Polar angle in degrees (0° = +z axis, 90° = xy-plane, 180° = -z axis)
    phi : float
        Azimuthal angle in degrees (0° = +x axis, 90° = +y axis)

    Returns:
    --------
    array : Field vector [x, y, z] as a unit vector

    Spherical to Cartesian conversion:
        x = sin(θ) cos(φ)
        y = sin(θ) sin(φ)
        z = cos(θ)
    """
    # Convert degrees to radians
    theta_rad = np.radians(theta)
    phi_rad = np.radians(phi)

    # Compute Cartesian components
    x = np.sin(theta_rad) * np.cos(phi_rad)
    y = np.sin(theta_rad) * np.sin(phi_rad)
    z = np.cos(theta_rad)

    return np.array([x, y, z])

# Example: Generate field vector
theta = 45  # 45° from z-axis
phi = 45.0    # 30° from x-axis in xy-plane

# compute field vector and print it
field_vector = generate_field_vector_from_theta_and_phi(theta, phi) * 0.1  # scale by field strength
print(f"Field vector for θ={theta}°, φ={phi}°: {field_vector}")

# ----------------------------
# Molecular geometry (ortho)
# ----------------------------
nitro_string = """
0 1
         C           -1.885946870148     1.189583649403    -0.119726128149
         C           -0.498436150781     1.207756989720    -0.095587674731
         C            0.177022900609     0.001244416404     0.003170781601
         C           -2.570320810975    -0.018930814056    -0.046076352491
         H           -2.433490753919     2.122014686561    -0.196502462215
         H            0.061742696521     2.131508866859    -0.151149504465
         H           -3.654760225882    -0.026965279238    -0.065539508188 
         N            1.653234989788     0.012187306043     0.029557107619  
         O            2.221533450190    -1.057060327273     0.115433821878  
         O            2.208198822578     1.089819744626    -0.036210769804    
         C           -1.871287950510    -1.217279026698     0.052462992898    
         H           -2.407539508613    -2.157651159151     0.109753013519    
         C           -0.483684066407    -1.215092211065     0.078168902711    
         H            0.087632738087    -2.130465296765     0.154596572701     
units angstrom
no_reorient
no_com
symmetry c1
"""


# ----------------------------
# Psi4 options
# ----------------------------
psi4_options = {
    "basis": "6-311G*",
    "scf_type": "df",          # density fitting
    "e_convergence": 1e-12,
    "d_convergence": 1e-12,
}

psi4.set_memory("24 GB")
psi4.core.set_output_file("nitrobenzene.out", False)
omega = 0.06615  # cavity frequency in atomic units (corresponding to ~1.8 eV)

# ----------------------------
# Build calculator
# ----------------------------
calculator = CQEDRHFCalculator(
    lambda_vector=field_vector, # molecule_string=ortho_string, #<-- CQEDRHFCalculator doesn't take molecule_string currently?
    psi4_options=psi4_options,
    omega=omega,
    density_fitting=True,
    charge=0,
    multiplicity=1
)

# ----------------------------
# Build orientation tracker
# ----------------------------
# We need initial coords + symbols for setup
mol = psi4.geometry(nitro_string)
psi4.set_options(psi4_options)
# compute psi4 rhf energy for initial geometry
e_rhf_1 = psi4.energy("scf", molecule=mol)
print(f"Initial RHF energy: {e_rhf_1:.6f} Hartree")
# use calculator to run a single point and get initial energy and forces
e_cqed_1 = calculator.energy(nitro_string)
print(f"Initial CQED-RHF energy: {e_cqed_1:.6f} Hartree")


# loop over theta from 0 to 180 and phi from 0 to 360 and compute energy for each orientation in increments of 2 degree
# print theta, phi, x, y, z, energy to a file or stdout in a formatted way
theta_list = np.arange(0, 181, num_theta_vals)
phi_list = np.arange(0, 361, num_phi_vals)


# Note that when theta = 0, all phi values give the same field vector (0, 0, 0.1) since the field is pointing along the z-axis. Similarly, when theta = 180, all phi values give the same field vector (0, 0, -0.1). The variation in energy will be more pronounced at intermediate theta values where the field vector changes direction in the xy-plane as phi varies.
# so we don't need to repeat these calculation for all phi values at theta = 0 and theta = 180, we can just compute it once for each of these theta values and print the result. This will save computational time while still capturing the key orientations of the field.
field_vector_theta_0 = generate_field_vector_from_theta_and_phi(0, 0) * 0.1
field_vector_theta_180 = generate_field_vector_from_theta_and_phi(180, 0) * 0.1

# update calculator with field vector for theta = 0 and compute energy
calculator.lambda_vector = field_vector_theta_0
c_cqed_theta_0 = calculator.energy(nitro_string)

# update calculator with field vector for theta = 180 and compute energy
calculator.lambda_vector = field_vector_theta_180
c_cqed_theta_180 = calculator.energy(nitro_string)


# open file for writing
with open("nitrobenzene_field_scan_results.txt", "w") as f:
    for theta in theta_list:
        for phi in phi_list:
            field_vector = generate_field_vector_from_theta_and_phi(theta, phi) * 0.1  # scale by field strength

            if theta == 0:
                field_vector = field_vector_theta_0
                e_cqed = c_cqed_theta_0
            elif theta == 180:
                field_vector = field_vector_theta_180
                e_cqed = c_cqed_theta_180
            else:   
                calculator.lambda_vector = field_vector
                e_cqed = calculator.energy(nitro_string)
            print(f"{theta:3f} {phi:3f} {field_vector[0]: .4f} {field_vector[1]: .4f} {field_vector[2]: .4f} {e_cqed:.12f}")
            f.write(f"{theta:3f} {phi:3f} {field_vector[0]: .4f} {field_vector[1]: .4f} {field_vector[2]: .4f} {e_cqed:.12f}\n")

f.close()


