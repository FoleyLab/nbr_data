import os
import matplotlib.pyplot as plt
import matplotlib as mpl

# ============================================================
# USER-SELECTABLE PARAMETERS
# ============================================================

save_path = "qed_ccsd_relaxed_vs_unrelaxed.png"
dpi = 150
figsize = (12, 5)

# ============================================================
# MATPLOTLIB RC PARAMETERS
# ============================================================

mpl.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          14,
    'axes.labelsize':     15,
    'axes.titlesize':     16,
    'xtick.labelsize':    13,
    'ytick.labelsize':    13,
    'legend.fontsize':    13,
    'lines.linewidth':    2.0,
    'lines.markersize':   8,
    'axes.linewidth':     1.2,
    'xtick.major.width':  1.0,
    'ytick.major.width':  1.0,
    'figure.dpi':         100,
    'savefig.dpi':        dpi,
    'savefig.bbox':       'tight',
})

# ============================================================
# DATA: QED-CCSD(2,2)
# ============================================================

hartree_to_kcal = 627.509

# --- dir_70_31: ortho - meta ---
lam_70 = [0.02, 0.04, 0.06, 0.08, 0.10]

E70_meta_unrel = [-3007.950230257695239, -3007.929917952462802,
                  -3007.897291260979273, -3007.853787768177426,
                  -3007.800881464791473]
E70_ortho_unrel = [-3007.948429049345577, -3007.929196702133595,
                   -3007.898152059528911, -3007.856519246041444,
                   -3007.805605599262435]

E70_meta_rel = [-3007.953275660864165, -3007.933025244492455,
                -3007.900494831265405, -3007.857089575379177,
                -3007.804261923117338]
E70_ortho_rel = [-3007.951921044827031, -3007.932945529768858,
                 -3007.902268823728718, -3007.861072299290754,
                 -3007.810595004073093]

d70_unrel = [(o - m) * hartree_to_kcal for o, m in zip(E70_ortho_unrel, E70_meta_unrel)]
d70_rel   = [(o - m) * hartree_to_kcal for o, m in zip(E70_ortho_rel,   E70_meta_rel)]

# --- dir_65_78: para - meta ---
lam_65 = [0.02, 0.04, 0.06, 0.08, 0.10]

E65_meta_unrel = [-3007.950250950728332, -3007.929991439987589,
                  -3007.897432239883074, -3007.853995023260723,
                  -3007.801149191801869]
E65_para_unrel = [-3007.942563718139354, -3007.923269314159825,
                  -3007.892092671240789, -3007.850231448629529,
                  -3007.798976102031247]

lam_65_rel = [0.02, 0.04, 0.08, 0.10]

E65_meta_rel = [-3007.953302390458703, -3007.933168412878331,
                -3007.857833311406466, -3007.805791594904349]
E65_para_rel = [-3007.950309055921025, -3007.931957315436193,
                -3007.862181624413097, -3007.813133009622106]

d65_unrel = [(p - m) * hartree_to_kcal for p, m in zip(E65_para_unrel, E65_meta_unrel)]
d65_rel   = [(p - m) * hartree_to_kcal for p, m in zip(E65_para_rel,   E65_meta_rel)]

# ============================================================
# PLOT
# ============================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

# -- Panel 1: dir_70_31 (ortho - meta) --
ax1.plot(lam_70, d70_unrel, 'o-',
         color='#1f77b4', markerfacecolor='#1f77b4',
         markersize=8, label='Unrelaxed geometry')
ax1.plot(lam_70, d70_rel, 's--',
         color='#d62728', markerfacecolor='#d62728',
         markersize=8, label='Relaxed geometry')
ax1.axhline(y=0, color='gray', linewidth=0.8, linestyle=':')
ax1.set_xlabel('λ (a.u.)')
ax1.set_ylabel('ΔE (kcal/mol)')
ax1.set_title('Ortho − Meta  (θ=70°, φ=31°)')
ax1.legend()
ax1.grid(True, alpha=0.25)

# -- Panel 2: dir_65_78 (para - meta) --
ax2.plot(lam_65, d65_unrel, 'o-',
         color='#1f77b4', markerfacecolor='#1f77b4',
         markersize=8, label='Unrelaxed geometry')
ax2.plot(lam_65_rel, d65_rel, 's--',
         color='#d62728', markerfacecolor='#d62728',
         markersize=8, label='Relaxed geometry')
ax2.axhline(y=0, color='gray', linewidth=0.8, linestyle=':')
ax2.set_xlabel('λ (a.u.)')
ax2.set_ylabel('ΔE (kcal/mol)')
ax2.set_title('Para − Meta  (θ=65°, φ=78°)')
ax2.legend()
ax2.grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), save_path))
print(f"Saved {save_path}")
