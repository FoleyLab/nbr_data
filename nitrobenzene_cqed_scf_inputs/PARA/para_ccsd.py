import numpy as np
import psi4


# ----------------------------
# Molecular geometry (meta)
# ----------------------------
meta_string = """
1 1
C         -0.511618296797     1.244386024531     0.732140048697
C          0.856500593203     1.251903714531     0.717948218697
C          1.524118723203     0.024661924531     0.713927788697
H         -1.071804396797     2.172682314531     0.745925708697
H          1.436128963203     2.163921874531     0.712099008697
N          3.008539583203     0.046097104531     0.698823798697
O          3.575097303203    -1.082768165469     0.699174708697
O          3.542114363203     1.190870854531     0.689202018697
C         -0.475464946797    -1.253402765469     0.742118638697
H         -1.008574426797    -2.197377955469     0.762945788697
C          0.892227703203    -1.221407805469     0.728065818697
H          1.498048423203    -2.116244695469     0.729653208697
C         -1.267906576797    -0.015841805469     0.712127418697
H         -2.116520796797    -0.025293215469     1.403161498697
Br        -2.114966986797    -0.034666925469    -1.121874081303
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
psi4.core.set_output_file("psi4_cc.out", False)

# Build orientation tracker
# ----------------------------
# We need initial coords + symbols for setup
mol = psi4.geometry(meta_string)
psi4.set_options(psi4_options)
# compute psi4 rhf energy for initial geometry
e_rhf = psi4.energy("scf", moleucule=mol)
#e_rhf_1 = psi4.energy("ccsd", molecule=mol)
#print(f"CCSD energy: {e_rhf_1:.6f} Hartree")




