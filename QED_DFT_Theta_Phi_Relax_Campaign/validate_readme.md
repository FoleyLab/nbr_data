Done. I added the generalized validator here:
[validate_qed_dft_campaign.py](/Users/jfoley19/Code/nbr_data/QED_DFT_Theta_Phi_Relax_Campaign/validate_qed_dft_campaign.py)
It now runs from the campaign top-level, discovers folders like para_th18_ph108, meta_th..., and ortho_th..., checks the expected completed layout, validates cell.json, opt_status.json, optimized.xyz, and *_qed_ccsd_input.json, and only imports/runs Psi4/CQED for selected expensive recomputes.
Useful invocations:
python validate_qed_dft_campaign.py --skip-incomplete
python validate_qed_dft_campaign.py --recompute-cell para_th0_ph0
python validate_qed_dft_campaign.py --recompute-count 5 --seed 123
python validate_qed_dft_campaign.py --recompute-probability 0.02 --skip-incomplete
Fast validation result on the current para campaign:
