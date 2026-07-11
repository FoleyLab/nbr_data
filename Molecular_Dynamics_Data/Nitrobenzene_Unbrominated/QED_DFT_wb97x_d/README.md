# Nitrobenzene Cavity-Orientation Bromination Regioselectivity — Analysis Workflow

This directory analyzes how coupling to an optical cavity's vacuum field changes the
regioselectivity of electrophilic bromination of nitrobenzene, by comparing an ab initio
molecular dynamics (AIMD) trajectory of the free molecule against orientation-dependent
energy-difference surfaces for the three possible Wheland (arenium) intermediates: ortho,
para, and meta.

The core idea: the molecule's instantaneous orientation relative to the cavity polarization
is described by two angles, theta (polar) and phi (azimuthal). A separate, higher-level
calculation (originally QED-CCSD, currently pQED(50,10) with EOM-CCSD, see "Swapping in new intermediate-energy data" below) gives
the ground-state energy of each Wheland intermediate as a function of (theta, phi). As the
AIMD trajectory tumbles the molecule through orientation space, this workflow looks up/interpolates
those intermediate energies at each trajectory frame's (theta, phi) and asks which intermediate —
ortho, para, or meta — is most stabilized at that moment.

## What's in this directory

**Trajectory input data** (already present, not modified by these scripts):
- `nitrobenzene_direction_A_wb97x_d_4000_ts.xyz`, `..._B_...`, `..._C_...`, `..._D_...` — four
  AIMD trajectories (QED-DFT/wB97X-D forces, NVE, ~2.5 ps, 4000 frames each), one per starting
  orientation ("direction"). Each is an XYZ trajectory where every frame is preceded by a header
  line: `Step <n>  E=<Hartree>  phi=<deg>  theta=<deg>`.

**Intermediate energy-surface data** (the (theta, phi) grid of Wheland-intermediate energies):
- `isomer_Nel_49_Nph_10_total_energies.dat` — 90x90 (theta, phi) grid (Nel=50 electronic states (gs + 49 excited states),
  Nph=10 photon Fock states in the Pauli-Fierz Hamiltonian, per the `EOM_CCSD` project's
  energy-scan script). **This is what `ENERGY_FILE` is currently set to in every script**

- `isomer_Nel_49_Nph_10_total_energies_CS.dat` — a second 90x90 grid using CS transformed PF Hamiltonian, 
   nearly numerically identical to the file above (differences ~1e-8 Hartree) 

- `bad_QED_CCSD_22_Results.txt` — a 24x24 (theta, phi) grid. The `bad_` prefix denotes
   this is the QED-CCSD(2,2) data that had a bug; **do not use without checking with Jay first**.


- All three files share the column format:
  `theta  phi  Ex  Ey  Ez  Para_E  Ortho_E  Meta_E` (angles in degrees, energies in Hartree).

**Analysis scripts** (each hardcoded to one specific trajectory "direction" — edit the
`MD_FILE`/`ENERGY_FILE` constants at the top of a script to point it at a different trajectory
or surface file):
1. `plot_para_meta_trajectory_overlay.py` (formerly `surfplot.py`) — trajectory-over-surface
   overlay (Para-Meta only), direction D. Produces `traj_overlay_direction_D.png`.
2. `plot_ortho_para_trajectory_overlays.py` (formerly `surfplot_para_meta.py`) —
   trajectory-over-surface overlay, BOTH Ortho-Meta and Para-Meta, direction C. Produces
   `traj_overlay_ortho_meta_direction_C.png` and `traj_overlay_para_meta_direction_C.png`.
3. `plot_para_meta_timeseries.py` (formerly `plt.py`) — DeltaE_para-meta(t) + relative energy
   time series, bilinear interpolation, direction B. Produces `deltaE_vs_time_direction_B.png`.
4. `plot_ortho_para_timeseries.py` (formerly `plt_para_meta.py`) — DeltaE_ortho-meta(t) AND
   DeltaE_para-meta(t) + relative energy time series, direction C. Produces
   `deltaE_vs_time_direction_C.png`.
5. `plot_timeseries_with_dwell_times.py` (formerly `global_timeseries_plt.py`) — the most
   complete time-series script: same plot as `plot_ortho_para_timeseries.py` plus basin
   classification and dwell-time bookkeeping (prints how long the trajectory spent in each
   stabilized basin), direction D. Produces `deltaE_vs_time_direction_D.png`.

## Environment / dependencies

Python 3 with `numpy`, `pandas`, `matplotlib`, and `scipy` (`scipy.interpolate.RegularGridInterpolator`,
used by all five scripts to map trajectory (theta, phi) points onto the energy surface via
bilinear interpolation). No other project-specific packages are required. Install with, e.g.:
```
pip install numpy pandas matplotlib scipy
```

## How to run

Each script is standalone (no shared imports between them) and is run directly:
```
python plot_para_meta_trajectory_overlay.py
python plot_ortho_para_trajectory_overlays.py
python plot_para_meta_timeseries.py
python plot_ortho_para_timeseries.py
python plot_timeseries_with_dwell_times.py
```
There is no enforced order between them — each one independently reads one trajectory `.xyz`
file and one energy-surface file and writes its own PNG(s) (and, for
`plot_timeseries_with_dwell_times.py`, prints dwell times to the console).

To analyze a different trajectory direction, or a different intermediate-energy surface, edit
the `MD_FILE` / `ENERGY_FILE` string constants near the top of the relevant script(s) — there are
no command-line arguments currently.

## Swapping in new (different level-of-theory) intermediate-energy data

The intermediate-energy surface has just been recomputed at a different level of theory than
the original QED-CCSD calculation. To point the workflow at the new data:

1. Confirm the new data file has the same column layout as the existing surface files:
   `theta  phi  Ex  Ey  Ez  Para_E  Ortho_E  Meta_E`, angles in degrees, energies in Hartree,
   with a one-line header followed by a divider line (the scripts do
   `pd.read_csv(ENERGY_FILE, sep='\s+', skiprows=[1])`, so row 2 must be a non-data divider row,
   not real data).
2. Change the `ENERGY_FILE = "..."` constant at the top of each script you plan to run to the new
   file's name/path (currently `"isomer_Nel_49_Nph_10_total_energies.dat"` in all five scripts).
   This is the only required change — none of the downstream interpolation, plotting, or
   classification logic depends on the level of theory.
3. Re-generate the surface plots and re-run the time series/dwell-time scripts; compare the new
   figures and dwell times against the ones already in this directory as a sanity check (see next
   section).
4. If the new grid's (theta, phi) sampling density or bounds differ from the current 90x90 grid
   spanning theta in [0, 180] and phi in [0, 360], double check
   `plot_para_meta_trajectory_overlay.py` / `plot_ortho_para_trajectory_overlays.py`'s
   `.reshape(num_t, num_p)` calls still make sense (they assume the input rows are pre-sorted as
   a perfect rectangular grid; `plot_para_meta_timeseries.py` / `plot_ortho_para_timeseries.py` /
   `plot_timeseries_with_dwell_times.py` use `.pivot(...)` instead, which is robust to row order).

## How to critically evaluate this analysis

Before trusting the output and moving on, check the following:

**Trajectory validity (NVE energy conservation).** Since the AIMD run is NVE starting from zero
kinetic energy, total energy should stay flat over the ~2.5 ps trajectory except for a small
initial-equilibration transient (the scripts already skip the first 10 frames for this reason).
Plot the raw `e_md` column (or look at the "Nitrobenzene Rel. E" curve already drawn in the time
series figures) and check for drift or a runaway trend — that would indicate an integrator/force
problem upstream of this analysis, not a modeling artifact.

**(theta, phi) mapping and interpolation onto the surface grid.** Check:
- Whether the trajectory's raw theta values ever fall outside the surface grid's sampled range
  (currently 0-180 for the grids in this directory). All five scripts use
  `scipy.interpolate.RegularGridInterpolator` with `fill_value=None` to map trajectory (theta,
  phi) points onto the surface, which *extrapolates* linearly for any out-of-bounds point rather
  than raising an error — silently extrapolated points could look reasonable in a plot but be
  numerically meaningless.
- The theta "branch correction" applied in every script's `parse_md_data`
  (`theta -> 180 - theta` whenever the raw parsed theta exceeds 100 degrees). This is an
  intentional, empirical fix: theta/phi are computed upstream (outside this directory) as
  `np.degrees(np.arccos(cos_theta))`, and arccos's principal range ([0, 180]) combined with a
  sign/orientation ambiguity in the reference axis (and arccos's numerical sensitivity near
  cos_theta = +/-1) can produce discontinuous jumps in the raw parsed trajectory. The `>100`
  fold restores continuity empirically; it is not a claim about (theta, phi) surface symmetry.
  Still worth double-checking against the new dataset if the underlying angle-extraction code
  ever changes.
- All five scripts now use the same bilinear interpolation method (consolidated 2026-07-10;
  `plot_para_meta_trajectory_overlay.py` previously used a KDTree nearest-neighbor lookup for an
  energy column it doesn't actually plot), so results are consistent across scripts and there's
  no separate nearest-neighbor mapping to reconcile.

**Classification thresholds and dwell-time logic.** The scripts use a +/-5 kcal/mol threshold
(matching the project's stated classification scheme) to decide whether an intermediate is
"stabilized" or "destabilized" relative to Meta. Also confirm the dwell-time bookkeeping in
`calculate_dwell_times` — it groups consecutive same-state frames via `itertools.groupby`, so a
single frame that briefly flips state and back will register as two additional (short) dwell
periods; decide whether that is the intended behavior for reporting "residence time."

**Units, end to end.** Angles: degrees throughout. Energies: Hartree in the raw data files,
converted to kcal/mol (`AU_TO_KCAL = 627.509`) only for plotting/classification — the underlying
`.dat`/`.txt` files should never be re-saved in kcal/mol. Time: each recorded MD step is 25
atomic units of time apart, so the step index is converted to femtoseconds via
`25 * 0.0241888 fs/au = 0.60472` fs/step (hardcoded as `6.04721e-16 * 1e15` in each script).

**Regression check against original results.** The original (QED-CCSD-level) analysis found
representative ortho-stabilizing (theta=70, phi=31) and para-stabilizing (theta=63, phi=63)
starting orientations, with roughly ~350 fs of ortho residence, ~200 fs of para residence, and
only ~50-150 fs of (transient) meta residence per trajectory. After swapping in the new
intermediate-energy data, re-run `plot_timeseries_with_dwell_times.py` for each direction and
compare its printed dwell times and the qualitative shape of the DeltaE-vs-time figure against
these numbers. Some shift is expected and legitimate given the different level of theory — the
sanity check is whether the *qualitative* picture (which basin dominates, roughly how transient
meta is) still holds, not exact numerical agreement.
