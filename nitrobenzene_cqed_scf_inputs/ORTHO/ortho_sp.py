import psi4
import numpy as np
from cqed_scf.calculator import CQEDCalculator

ortho_string = """
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
1 1
units angstrom
no_reorient
no_com
symmetry c1
"""


# -----------------------------------------------------------------------------
# QED-DFT settings
# -----------------------------------------------------------------------------

PSI4_OPTIONS = {
    "basis": "6-31G",
    "reference": "rks",
    "scf_type": "df",
    "e_convergence": 1e-9,
    "d_convergence": 1e-9,
    "dft_radial_points": 99,
    "dft_spherical_points": 590,
    "dft_pruning_scheme": "none",
}

mol = psi4.geometry(ortho_string)
psi4.set_options(PSI4_OPTIONS)

e = psi4.energy("scf")


#LAMBDA_VECTOR = np.array([0.078, 0.055, 0.027])
#OMEGA = 0.06615
#FUNCTIONAL = "wb97x-d"

#----------------------------------------------------
# Create cqed-scf calculator with the settings above
#----------------------------------------------------
#calc = CQEDCalculator(
#    lambda_vector = LAMBDA_VECTOR,
#    psi4_options = PSI4_OPTIONS,
#    omega = OMEGA,
#    charge = 0,
#    multiplicity = 1,
#    density_fitting = True,
#    functional = FUNCTIONAL,
#    debug = False
#)

#calc.energy(ortho_string)

