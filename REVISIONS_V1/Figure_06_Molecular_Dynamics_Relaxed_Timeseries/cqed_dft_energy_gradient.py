"""
CQED-DFT energy and nuclear gradient using CQEDConfig and CQEDCalculator.
"""

import numpy as np
import psi4

from cqed_scf import CQEDCalculator, CQEDConfig


# ---------------------------------------------------------
# Geometry
# ---------------------------------------------------------

geometry = """
0 1
         C           -1.901922841579     1.172760393927    -0.047772892327
         C           -0.476703461579     1.208889463927    -0.052552462327
         C            0.199130338421     0.032560163927    -0.040813672327
         C           -2.617281111579    -0.030143496073    -0.014796322327
         H           -2.441941791579     2.113963523927    -0.059190102327
         H            0.061569388421     2.148636453927    -0.051942432327
         H           -3.698948961579    -0.020371166073    -0.004154112327
         N            1.670887158421     0.035666643927     0.020806147673
         O            2.218059368421    -1.081966866073     0.128768557673
         O            2.243651238421     1.138122873927    -0.039099092327
         C           -1.926962411579    -1.227400606073     0.005173817673
         H           -2.453758301579    -2.175200806073     0.038655707673
         C           -0.458335831579    -1.276806306073    -0.036936782327
         H            0.022343538421    -2.023338446073     0.599315217673
no_reorient
no_com
symmetry c1
"""


# ---------------------------------------------------------
# Psi4 options
# ---------------------------------------------------------

psi4.set_memory("4 GB")

psi4_options = {
    "basis": "6-311g*",
    "scf_type": "df",
    "e_convergence": 1e-10,
    "d_convergence": 1e-9,
    "dft_radial_points": 99,
    "dft_spherical_points": 590,
    "dft_pruning_scheme": "none",
}


# ---------------------------------------------------------
# Build CQED configuration
# ---------------------------------------------------------
field_vector = np.array([0.07878123598, 0.0551632153, 0.02739592187])
omega = 0.06615

config_cav = CQEDConfig(
    lambda_vector=field_vector,
    omega=omega,
    psi4_options=psi4_options,
    reference="rks",
    functional="wb97x-d",
    density_fitting=True,
    charge=0,
    multiplicity=1,
    dispersion_policy="post_scf",
    debug=False,
    quiet=False,  # NORMAL: default verbosity.
    # NOTE: gradient-path prints (CQEDGradient timings/components) are not yet
    # routed through cqed_scf.output (Stage B), so quiet=True would silence only
    # the SCF/energy portion. Normal mode is the meaningful setting here today.
)

config_no_cav = CQEDConfig(
    lambda_vector=np.array([0.0, 0.0, 0.0]),
    omega=omega,
    psi4_options=psi4_options,
    reference="rks",
    functional="wb97x-d",
    density_fitting=True,
    charge=0,
    multiplicity=1,
    dispersion_policy="post_scf",
    debug=False,
    quiet=False,  # NORMAL: default verbosity.
    # NOTE: gradient-path prints (CQEDGradient timings/components) are not yet
    # routed through cqed_scf.output (Stage B), so quiet=True would silence only
    # the SCF/energy portion. Normal mode is the meaningful setting here today.
)


# ---------------------------------------------------------
# Run CQED-DFT energy + gradient
# ---------------------------------------------------------

calc_cav = CQEDCalculator(config=config_cav)
calc_no_cav = CQEDCalculator(config=config_no_cav)
energy_cav  = calc_cav.energy(geometry)
energy_no_cav = calc_no_cav.energy(geometry)


print("\nCQED-DFT energy and gradient")
print("============================")
print(f"Reference   : {config_cav.reference}")
print(f"Functional  : {config_cav.functional}")
print(f"Energy w/ Cavity      : {energy_cav:.12f} Eh")
print(f"Energy w/o Cavity     : {energy_no_cav:.12f} Eh")
print(f"Energy Difference (kcal/mol)     : {(energy_cav - energy_no_cav) * 627.509:.12f} kcal/mol")
