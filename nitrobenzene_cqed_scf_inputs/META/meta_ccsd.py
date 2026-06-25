import numpy as np
import psi4


# ----------------------------
# Molecular geometry (meta)
# ----------------------------
meta_string = """
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
e_rhf_1 = psi4.energy("ccsd", molecule=mol)
print(f"CCSD energy: {e_rhf_1:.6f} Hartree")




