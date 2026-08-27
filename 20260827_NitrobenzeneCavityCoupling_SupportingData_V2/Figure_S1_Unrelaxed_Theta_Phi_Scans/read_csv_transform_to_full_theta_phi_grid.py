import numpy as np

# --- Configuration ---
#INPUT_FILE = "intermediate_scans.csv"
INPUT_FILE = "unrelaxed_qed_ccsd_intermediate_scans.csv"
CCSD_OUTPUT = "unrelaxed_QED_CCSD_theta_phi_scan.txt"
SCF_OUTPUT = "unrelaxed_QED_SCF_theta_phi_scan.txt"
# Exact headers and divider requested
HEADER_1 = "     theta        phi           Ex           Ey           Ez               Para_E              Ortho_E               Meta_E"
HEADER_2 = "---------------------------------------------------------------------------------------------------------------------------"

def generate_field_vector_from_theta_and_phi(theta, phi):
    """Generate a unit field vector from spherical coordinate transformation."""
    theta_rad = np.radians(theta)
    phi_rad = np.radians(phi)
    
    x = np.sin(theta_rad) * np.cos(phi_rad)
    y = np.sin(theta_rad) * np.sin(phi_rad)
    z = np.cos(theta_rad)
    
    return np.array([x, y, z])

def build_symmetry_lookups():
    # Load raw data from CSV
    data = np.genfromtxt(INPUT_FILE, skip_header=1, delimiter=",")
    theta_raw, phi_raw = data[:, 0], data[:, 1]
    
    e_scf_ortho = data[:, 2]
    e_ccsd_ortho = data[:, 3]
    e_scf_meta = data[:, 4]
    e_ccsd_meta = data[:, 5]
    e_scf_para = data[:, 6]
    e_ccsd_para = data[:, 7]
    
    # Create dictionaries for fast coordinate lookup
    ccsd_lookup = {}
    scf_lookup = {}
    
    # Track anchor values for polar degeneracies
    # We will grab the first available row matching theta=0 or theta=180
    ccsd_theta0, scf_theta0 = None, None
    ccsd_theta180, scf_theta180 = None, None
    
    for i in range(len(data)):
        t_val = theta_raw[i]
        p_val = phi_raw[i]
        t_key = round(t_val, 2)
        p_key = round(p_val, 2)
        
        energies_ccsd = (e_ccsd_para[i], e_ccsd_ortho[i], e_ccsd_meta[i])
        energies_scf = (e_scf_para[i], e_scf_ortho[i], e_scf_meta[i])
        
        ccsd_lookup[(t_key, p_key)] = energies_ccsd
        scf_lookup[(t_key, p_key)] = energies_scf
        
        if np.isclose(t_val, 0.0) and ccsd_theta0 is None:
            ccsd_theta0 = energies_ccsd
            scf_theta0 = energies_scf
        if np.isclose(t_val, 180.0) and ccsd_theta180 is None:
            ccsd_theta180 = energies_ccsd
            scf_theta180 = energies_scf

    # FIX: when the input scan only covers the canonical theta in [0,90] half (as
    # intermediate_scans.csv does), there is no raw theta=180 row, so ccsd_theta180/
    # scf_theta180 stay None here. Previously that meant fetch_energy_with_symmetry's
    # Rule 2 never fired for theta=180, and nearly every phi at that pole silently fell
    # through to the (0.0, 0.0, 0.0) placeholder -- only phi=180 accidentally worked,
    # because that's the one case where the Rule 5 inversion lookup happens to land on
    # the single stored theta=0 point.
    #
    # The theta=0 and theta=180 poles must carry the identical (phi-independent) energy:
    # the field vector (Ex,Ey,Ez) is phi-independent at both poles (sin(theta)=0 there),
    # and the antipodal symmetry E(theta,phi) = E(180-theta, phi+180 mod 360) maps the
    # whole theta=0 ring onto the whole theta=180 ring. So if no explicit theta=180 data
    # was scanned, derive it directly from the theta=0 anchor instead of leaving it None.
    if ccsd_theta180 is None:
        ccsd_theta180 = ccsd_theta0
    if scf_theta180 is None:
        scf_theta180 = scf_theta0

    # Dynamically extract the resolution step sizes
    d_theta = np.diff(np.unique(theta_raw)).min()
    d_phi = np.diff(np.unique(phi_raw)).min()

    return ccsd_lookup, scf_lookup, d_theta, d_phi, ccsd_theta0, scf_theta0, ccsd_theta180, scf_theta180

def fetch_energy_with_symmetry(t, p, lookup, theta0_data, theta180_data):
    """Applies contextual polar, azimuthal, and inversion symmetries to find energy."""
    # Rule 1: Polar Degeneracy at Theta = 0
    if np.isclose(t, 0.0) and theta0_data is not None:
        return theta0_data
        
    # Rule 2: Polar Degeneracy at Theta = 180
    if np.isclose(t, 180.0) and theta180_data is not None:
        return theta180_data
        
    # Rule 3: Azimuthal Wrap Symmetry (phi = 360 maps to phi = 0)
    if np.isclose(p, 360.0):
        p = 0.0
        
    t_key, p_key = round(t, 2), round(p, 2)
    
    # Rule 4: Match exact coordinate point if it exists
    if (t_key, p_key) in lookup:
        return lookup[(t_key, p_key)]
    
    # Rule 5: Vector Inversion Symmetry (-vec fallback)
    t_inv = round(180.0 - t, 2)
    p_inv = round((p + 180.0) % 360.0, 2)
    if p_inv == 360.0:
        p_inv = 0.0
        
    if (t_inv, p_inv) in lookup:
        return lookup[(t_inv, p_inv)]

    # If every rule above fails, the input scan has an actual coverage gap -- fail loudly
    # instead of silently writing a (0.0, 0.0, 0.0) placeholder that looks like real data.
    raise ValueError(
        f"No data or symmetry-derived value found for (theta={t}, phi={p}); "
        f"checked direct match at (theta={t_key}, phi={p_key}) and inversion match at "
        f"(theta={t_inv}, phi={p_inv}). The input scan may not fully cover the canonical "
        f"half-domain needed to reconstruct this point."
    )

def convert_data():
    # 1. Build lookup tables and extract boundary anchors
    (ccsd_lookup, scf_lookup, d_theta, d_phi, 
     ccsd_t0, scf_t0, ccsd_t180, scf_t180) = build_symmetry_lookups()
    
    # 2. Define the full desired coordinate loops
    full_theta = np.arange(0, 180 + d_theta/2, d_theta)
    full_phi = np.arange(0, 360 + d_phi/2, d_phi)
    
    # 3. Open files and write out structured fixed-width formatting
    with open(CCSD_OUTPUT, "w") as f_ccsd, open(SCF_OUTPUT, "w") as f_scf:
        f_ccsd.write(f"{HEADER_1}\n{HEADER_2}\n")
        f_scf.write(f"{HEADER_1}\n{HEADER_2}\n")
        
        for t in full_theta:
            for p in full_phi:
                # Generate directional unit vectors
                vec = generate_field_vector_from_theta_and_phi(t, p)
                ex, ey, ez = vec[0], vec[1], vec[2]
                
                # Fetch energies enforcing boundary and inversion symmetries
                para_ccsd, ortho_ccsd, meta_ccsd = fetch_energy_with_symmetry(t, p, ccsd_lookup, ccsd_t0, ccsd_t180)
                para_scf, ortho_scf, meta_scf = fetch_energy_with_symmetry(t, p, scf_lookup, scf_t0, scf_t180)
                
                # Write to text using the exact fixed-width layout structure
                line_ccsd = f"   {t:10.1f} {p:11.1f} {ex:12.7f} {ey:12.7f} {ez:12.7f}   {para_ccsd:21.14f}   {ortho_ccsd:21.14f}   {meta_ccsd:21.14f}\n"
                line_scf  = f"   {t:10.1f} {p:11.1f} {ex:12.7f} {ey:12.7f} {ez:12.7f}   {para_scf:21.14f}   {ortho_scf:21.14f}   {meta_scf:21.14f}\n"
                
                f_ccsd.write(line_ccsd)
                f_scf.write(line_scf)
                
    print(f"Successfully generated files with boundary degeneracies:\n - {CCSD_OUTPUT}\n - {SCF_OUTPUT}")

if __name__ == "__main__":
    convert_data()
