import json
import re

import numpy as np
import psi4
from cqed_scf.calculator import CQEDCalculator
from cqed_scf.utils import generate_field_vector_from_theta_and_phi

# ---------------------------------------------------------------------------
# 1. Load parameters from cell.json
# ---------------------------------------------------------------------------
with open("cell.json") as f:
    cell = json.load(f)

theta = cell["theta"]
phi = cell["phi"]
LAMBDA_MAG = cell["magnitude"]
OMEGA = 0.06615
FUNCTIONAL = "wb97x"

FIELD_VECTOR = generate_field_vector_from_theta_and_phi(theta=theta, phi=phi)
LAMBDA_VECTOR = LAMBDA_MAG * FIELD_VECTOR
print(f"Loaded parameters from cell.json: theta={theta}, phi={phi}, magnitude={LAMBDA_MAG}")
print(f"Computed FIELD_VECTOR: {FIELD_VECTOR}")
print(f"Computed LAMBDA_VECTOR: {LAMBDA_VECTOR}")

# Verify that FIELD_VECTOR * LAMBDA_MAG matches cell.json lambda_vector
expected_lv = np.array(cell["lambda_vector"])
computed_lv = np.array(LAMBDA_VECTOR)
np.testing.assert_allclose(
    computed_lv, expected_lv, atol=1e-15,
    err_msg="lambda_vector from cell.json does not match LAMBDA_MAG * FIELD_VECTOR"
)
print("✓ cell.json lambda_vector verified against theta/phi/magnitude")

# ---------------------------------------------------------------------------
# 2. Parse the optimized geometry and expected QA values from optimized.xyz
# ---------------------------------------------------------------------------
with open("optimized.xyz", "r") as f:
    xyz_lines = f.readlines()

n_atoms = int(xyz_lines[0].strip())
comment = xyz_lines[1].strip()

# Parse |g_proj| and E from the comment line
m = re.search(r'\|\s*g_proj\s*\|\s*=\s*([\d.eE+-]+)', comment)
_EXPECTED_G_PROJ_NORM = float(m.group(1)) if m else None
m = re.search(r'E\s*=\s*([\d.eE+-]+)', comment)
_EXPECTED_ENERGY = float(m.group(1)) if m else None

# Build geometry string for Psi4
geom_lines = []
for line in xyz_lines[2:2 + n_atoms]:
    parts = line.strip().split()
    if len(parts) >= 4:
        elem, x, y, z = parts[0], parts[1], parts[2], parts[3]
        geom_lines.append(f"{elem:4s} {x:>20s} {y:>20s} {z:>20s}")

geometry_string = f"""1 1
{chr(10).join(geom_lines)}
units angstrom
no_reorient
no_com
symmetry c1
"""
print(f"Parsed optimized geometry with {n_atoms} atoms from optimized.xyz")
print(f"Expected |g_proj|: {_EXPECTED_G_PROJ_NORM}, Expected Energy: {_EXPECTED_ENERGY}")
print(f"Geometry string for Psi4:\n{geometry_string}")


# ---------------------------------------------------------------------------
# 4. Run the DFT calculation and compare with expected values
# ---------------------------------------------------------------------------
PSI4_OPTIONS = {
    "basis": "6-311G*",
    "reference": "rks",
    "scf_type": "df",
    "e_convergence": 1e-9,
    "d_convergence": 1e-9,
}

calc = CQEDCalculator(
    lambda_vector=LAMBDA_VECTOR,
    psi4_options=PSI4_OPTIONS,
    omega=OMEGA,
    density_fitting=True,
    charge=1,
    multiplicity=1,
    functional=FUNCTIONAL,
    reference="rks",
    dispersion_policy="none",
    debug=False,
)

energy, gradient, _ = calc.energy_and_projected_gradient(geometry_string)
g_proj_norm = np.linalg.norm(gradient)

print(f"Computed energy:         {energy:.12f} Hartree")
print(f"Expected energy:         {_EXPECTED_ENERGY:.12f} Hartree")
print(f"Computed |g_proj|:       {g_proj_norm:.6e} Hartree/Bohr")
print(f"Expected |g_proj|:       {_EXPECTED_G_PROJ_NORM:.6e} Hartree/Bohr")

# Validate against expected values from optimized.xyz comment
assert abs(energy - _EXPECTED_ENERGY) < 1e-8, \
    f"Energy mismatch: computed {energy}, expected {_EXPECTED_ENERGY}"
assert abs(g_proj_norm - _EXPECTED_G_PROJ_NORM) < 1e-6, \
    f"|g_proj| mismatch: computed {g_proj_norm}, expected {_EXPECTED_G_PROJ_NORM}"

print("\n✓ All automated checks passed.")
