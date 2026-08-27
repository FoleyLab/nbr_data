# Grid campaign (no frequencies): theta/phi-grid constrained QED-DFT relaxations

This document covers `grid_campaign_no_freq.py`, `ccsd_input.py`, and the small
addition to `campaign.py` that support them -- what they do, how to run them,
and the reasoning behind how they're built. It assumes familiarity with the
existing `campaign.py` / `trajectory_campaign.py` / `cavity_common.py` trio;
where this campaign differs from those, it says so explicitly.

## What this campaign does

For one Wheland-intermediate isomer at a time (ortho, meta, or para), it runs
a constrained QED-DFT geometry relaxation at every (theta, phi) cavity
orientation listed in `intermediate_scans.csv` -- the same grid used for the
QED-CCSD/EOM-CCSD energy scan (`Energy_Scan_Data/QED_CCSD/intermediate_scans.csv`
in the data repo). Only the `theta,phi` columns of that CSV are read; its
energy columns are a different level of theory's precomputed scan and are not
consulted here.

Unlike `trajectory_campaign.py` / `trajectory_campaign_no_freq.py` -- where
every (theta, phi) direction is an independent cell, safe to run in any order
or fully in parallel -- this campaign **warm-starts across phi within a fixed
theta**: phi=0 is optimized starting from the isomer's pristine `.xyz`; the
next phi in that theta ring is then optimized starting from phi=0's *converged*
geometry; and so on through every phi listed for that theta. Each theta ring
starts fresh from the pristine isomer geometry at phi=0 -- it is never seeded
from a neighboring theta ring -- so the 11 theta rings (and the 3 isomers) are
all independent of each other and safe to run concurrently. Only the phi steps
*within* one ring are sequential.

There is no frequency/Hessian step at all: this is a "light relaxation," gated
on the projected-gradient norm alone (see Convergence below), not a full
optimize-then-vibrate pipeline. On every cell that converges, a
`<isomer>_th<theta>_ph<phi>_qed_ccsd_input.json` is also written, ready to
hand to a downstream QED-CCSD calculation at that same (theta, phi) direction.
It's named after the cell rather than just `qed_ccsd_input.json` on purpose:
these files are meant to be gathered up and run on a different machine, by a
different person, later -- a bare `qed_ccsd_input.json` would be ambiguous
the moment more than one of them sits in the same place (e.g. after copying
several off of `runs_grid_no_freq/` into one folder to launch), so the
isomer/theta/phi identity travels with the filename itself.

## Methodology: why a warm-start chain

The rationale is physical continuity, not just speed: at fixed theta, as phi
sweeps the cavity field around the z-axis, the field vector rotates smoothly,
and the optimal (lowest-energy, field-coupled) geometry should follow it
smoothly too. Starting each phi step from the *previous* phi's converged
geometry, rather than from the pristine isomer every time, gives BFGS a much
better initial guess and should reach the (loose) convergence gate in fewer
restart cycles per cell than 191 independent from-scratch optimizations would.

This only holds within one theta ring, which is why theta rings don't chain
into each other: two different theta values are not "adjacent" in the same
continuous-path sense that two consecutive phi values are (jumping from one
theta ring to the next isn't a small perturbation of the cavity direction in
the way one phi step is), so seeding across theta wasn't assumed to carry the
same benefit, and doing so wasn't asked for. Each ring instead re-anchors at
the isomer's own pristine geometry at phi=0, so a bad chain in one ring can't
propagate into another.

**Consequence for parallelism.** Because phi steps within a ring are strictly
sequential and depend on each other, the unit of parallel work is the
**(isomer, theta) ring**, not the individual cell -- this is the one place
this campaign's execution model differs structurally from
`trajectory_campaign.py`'s. `--list-groups` / `--only-group` expose rings
(rather than individual cell ids) for a job launcher, for exactly this reason.

**Consequence for failure handling.** If a phi step doesn't converge within
`MAX_RESTARTS` (imported unchanged from `campaign.py`), that ring **halts**:
later phi values in the same ring are left pending rather than warm-started
from a geometry that never actually converged. This was a deliberate choice
(over silently skipping the stalled cell and warm-starting the next one from
the last *good* geometry, or falling back to the pristine isomer) -- it keeps
every warm-start in the campaign provably seeded from a converged geometry,
at the cost of needing a human to look at (and likely re-run, possibly with
adjusted settings) a stalled cell before its ring can finish. Re-running the
campaign after a stall is fixed resumes the ring from wherever it left off
(via the same per-cell `DONE`-marker check described in Output layout below)
-- it does not restart the ring from phi=0.

## Convergence: a looser, separate gate

`campaign.py`'s existing `CONV_THRESHOLD` (2.5e-3) is shared by
`campaign.py` and `trajectory_campaign.py`'s frequency-bound cells. This
campaign wants a different, explicitly *looser-than-that-was-meant-to-be*
gate for a lighter relaxation: `GRID_CONV_THRESHOLD = 2e-3` in
`grid_campaign_no_freq.py`.

Rather than repoint the shared `CONV_THRESHOLD`, which would silently change
behavior for the existing frequency campaigns, `campaign.optimize_cell` grew
an optional `conv_threshold=CONV_THRESHOLD` parameter (default unchanged, so
`campaign.py` and `trajectory_campaign.py` call it exactly as before). This
campaign is the one caller that passes `conv_threshold=GRID_CONV_THRESHOLD`
explicitly. The value actually used is still recorded per cell in
`opt_status.json["conv_threshold"]`, same as before -- provenance doesn't
depend on remembering which campaign produced a given `runs*/` directory.

## Files

- `intermediate_scans.csv` -- a snapshot of the (theta, phi) scan grid, copied
  in from the data repo so this campaign is self-contained and reproducible
  even if the original copy changes later. Only `theta` and `phi` are used.
- `ccsd_input.py` (new) -- builds the QED-CCSD input JSON for one converged
  cell. Deliberately a separate module from `cavity_common.py`: the latter is
  the single source of truth for the QED-DFT *optimization* physics that this
  campaign (and `campaign.py` / `trajectory_campaign.py`) actually run;
  `ccsd_input.py` only formats the *output* of that optimization into an
  input file for a different downstream program (QED-CCSD), at a separate,
  fixed lambda-magnitude convention (0.1, independent of whatever magnitude
  the optimization itself used). Keeping the two apart means a change to the
  DFT opt options can never accidentally perturb the CCSD input template, and
  a change to the CCSD input format can never perturb the optimization.
  Reuses `cavity_common.lambda_vector_for(theta, phi, magnitude=1.0)` for the
  `qed_polvecs` unit direction, and `cavity_common.OMEGA` / `CHARGE` /
  `MULTIPLICITY` / `PSI4_OPTIONS["basis"]` for the fields the two programs
  must actually agree on -- so there is still exactly one place each of those
  physical values is defined, even across the module boundary.
  `cavity_common` is imported lazily, inside the functions that need it (not
  at module level), matching `campaign.py`'s own convention: importing
  `ccsd_input.py` to build a manifest or run a `--dry-run` shouldn't require
  psi4/ase/cqed_scf to be installed.
- `grid_campaign_no_freq.py` (new) -- the driver: parses the scan file,
  builds the per-isomer per-theta phi chains, runs them with the warm-start /
  halt-on-stall logic described above, and writes
  `<cell_id>_qed_ccsd_input.json` on every convergence.
- `campaign.py` -- unchanged except for `optimize_cell`'s new optional
  `conv_threshold` parameter.

## How to run

One isomer per invocation, same convention as `trajectory_campaign.py`:

```
python grid_campaign_no_freq.py --isomer ortho                 # run/resume every ring
python grid_campaign_no_freq.py --isomer ortho --dry-run        # print the manifest
python grid_campaign_no_freq.py --isomer ortho --list-groups    # ring ids, for a launcher
python grid_campaign_no_freq.py --isomer ortho --only-group th63   # just one ring
python grid_campaign_no_freq.py --isomer ortho --summarize      # rebuild the rollup CSV/JSON
```

Embarrassingly-parallel run, sharded by ring (NOT by individual cell, since
phi steps within a ring must stay sequential in one process). `--list-groups`
prints one ring id per line, so piping it through `parallel -j N` launches
one process per ring, up to N concurrently:

```
export PSI4_THREADS=4
python grid_campaign_no_freq.py --isomer ortho --list-groups | \
    parallel -j 4 python grid_campaign_no_freq.py --isomer ortho --only-group {} --no-summarize
python grid_campaign_no_freq.py --isomer ortho --summarize
```

`--no-summarize` on the per-ring shards avoids every parallel process racing
to rewrite the same rollup file; the final `--summarize` call rebuilds it
once, after every ring has finished.

### Worked example: ortho, all 11 theta rings at once, 11-way concurrency

`ortho` has exactly 11 theta rings (theta = 0, 9, 18, ..., 90), so if the
machine can run 11 jobs at once, every ring can go in parallel in a single
pass -- set `-j 11` (or `-j0`, GNU parallel's "no limit" shorthand, since
`--list-groups` for one isomer never emits more than 11 lines) and each of
the 11 `--only-group` invocations gets its own ring, all starting together:

```
cd opt_only_campaign
export PSI4_THREADS=4   # tune per-job thread count so 11 * PSI4_THREADS <= cores available

# 1. Sanity-check the manifest first (no compute, just prints it):
python grid_campaign_no_freq.py --isomer ortho --dry-run

# 2. See the 11 ring ids that --list-groups will feed to parallel:
python grid_campaign_no_freq.py --isomer ortho --list-groups
#   ortho_th0
#   ortho_th9
#   ortho_th18
#   ortho_th27
#   ortho_th36
#   ortho_th45
#   ortho_th54
#   ortho_th63
#   ortho_th72
#   ortho_th81
#   ortho_th90

# 3. Launch all 11 rings concurrently, one process each:
python grid_campaign_no_freq.py --isomer ortho --list-groups | \
    parallel -j 11 \
    python grid_campaign_no_freq.py --isomer ortho --only-group {} --no-summarize

# 4. Once all 11 have finished (or halted on a stall -- check the console
#    output / runs_grid_no_freq/*/DONE for which), rebuild the rollup once:
python grid_campaign_no_freq.py --isomer ortho --summarize
```

If a ring halts on a stall (see Methodology above), the other 10 are
unaffected and keep running to completion; re-running step 3 afterward only
re-touches the halted ring's still-pending cells (everything else is already
`DONE` and gets skipped), so it's safe to just re-run the same `parallel`
line again once the stall is investigated. Repeat the whole sequence with
`--isomer meta` and `--isomer para` for the other two intermediates.

## Output layout

Each cell lives in `runs_grid_no_freq/<isomer>_th<theta>_ph<phi>/`:

- `cell.json` -- the cell's manifest entry plus its derived lambda vector.
- `opt_traj.xyz` -- every BFGS restart attempt's geometry, appended.
- `opt_status.json` -- converged flag, attempt count, gradient-norm history,
  final energy, `conv_threshold` used, wall time.
- `optimized.xyz` -- final converged geometry (only written if converged;
  this is what seeds the next phi in the ring).
- `<isomer>_th<theta>_ph<phi>_qed_ccsd_input.json` -- e.g.
  `ortho_th63_ph126_qed_ccsd_input.json`; only written on convergence (see
  below). Named after the cell, not just `qed_ccsd_input.json`, so it stays
  identifiable once copied out of its `runs_grid_no_freq/<cell_id>/` folder.
- `DONE` -- marker file; its presence is what makes a cell skippable on
  re-run.

Per-isomer rollups: `runs_grid_no_freq/grid_campaign_no_freq_<isomer>_status.{csv,json}`,
one row per cell (`cell_id`, `group_id`, `theta`, `phi`, `chain_index`,
convergence/gradient/energy columns, and whether the CCSD input was written).

## QED-CCSD input generation

`<isomer>_th<theta>_ph<phi>_qed_ccsd_input.json` mirrors the example input
Jay provided:

- `geometry.coordinates` -- the converged geometry for this cell.
- `SCF.qed_polvecs` -- the unit field-direction vector for this cell's
  (theta, phi), from `generate_field_vector_from_theta_and_phi` via
  `cavity_common.lambda_vector_for(theta, phi, magnitude=1.0)`.
- `SCF.qed_lambdas` -- always `[0.1]`, the fixed coupling magnitude for every
  cell in this campaign (independent of the DFT optimization's own lambda
  magnitude convention, should they ever differ).
- `SCF.qed_omegas`, `SCF.charge`, `SCF.multiplicity`, `basis.basisset` --
  pulled from `cavity_common` rather than re-hardcoded, so they can't drift
  out of sync with whatever the optimization actually used.
- Everything else (`common`, the rest of `SCF`, `CD`, `CC`, `TASK`,
  `geometry.units`/`noorient`) is fixed, taken literally from the example.

## Verification

psi4/ase/cqed_scf aren't available in the environment this was developed in,
so the actual DFT optimization couldn't be run end-to-end here. What *was*
verified directly, with `campaign.optimize_cell` and `cavity_common` mocked
out (a fake optimizer that deterministically "converges" or "stalls" cells on
demand, and a fake `cavity_common` module supplying the same lambda-vector
math): manifest/ring construction against the real `intermediate_scans.csv`
(11 rings, 191 cells, matching the CSV's actual theta/phi counts exactly,
including the pole's single phi and the equator ring's half-range phi);
warm-start seeding (numerically confirmed a later phi's geometry is derived
from the *previous* phi's converged geometry, not the pristine isomer);
halt-on-stall (a forced stall partway through a ring leaves later phi
untouched); resume-after-fixing-a-stall (the ring picks up and completes,
without re-running the already-converged earlier phi); idempotent re-runs;
the per-cell-named `<cell_id>_qed_ccsd_input.json` structure and unit-norm
`qed_polvecs`; and the `summarize()` rollup. Before running the real campaign, it's still worth
doing a `--dry-run` (to sanity-check the manifest against the current
`intermediate_scans.csv`) and a small `--only-group` smoke test on hardware
with psi4/cqed_scf installed.

## Relationship to the existing campaigns

- `campaign.py` -- unchanged in behavior; only gained the optional
  `conv_threshold` parameter on `optimize_cell` that this campaign uses.
- `trajectory_campaign.py` / `trajectory_campaign_no_freq.py` -- the pattern
  this campaign's per-cell bookkeeping (DONE markers, `opt_status.json`,
  rollup CSV/JSON) follows closely. The structural difference is the
  execution model: those treat every direction as independent and
  parallelize per-cell; this campaign chains phi within a theta ring and
  parallelizes per-ring.
