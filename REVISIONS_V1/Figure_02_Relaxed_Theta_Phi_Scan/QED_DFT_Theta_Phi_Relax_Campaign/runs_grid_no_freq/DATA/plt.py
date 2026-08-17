import pandas as pd
import numpy as np

# read grid_campaign_no_freq_opt_energies_and_gnorm.csv into data frame
df = pd.read_csv('grid_campaign_no_freq_opt_energies_and_gnorm.csv')

# print the keys
print(df.keys())

# keys are
#Index(['theta', 'phi', 'E_opt_ortho', 'E_opt_meta', 'E_opt_para',
#       'g_norm_ortho', 'g_norm_meta', 'g_norm_para'],
#      dtype='str')

# find maximum value of E_opt_ortho - E_opt_meta and its corresponding theta and phi
df['E_diff'] = (df['E_opt_ortho'] - df['E_opt_meta']) * 627.5094740631  # convert hartree to kcal/mol
max_E_diff = df['E_diff'].min()
max_E_diff_row = df[df['E_diff'] == max_E_diff]
max_theta = max_E_diff_row['theta'].values[0]
max_phi = max_E_diff_row['phi'].values[0]

# find maximum value of g_norm_ortho - g_norm_meta 
df['g_norm_diff'] = df['g_norm_ortho'] - df['g_norm_meta']
max_g_norm_diff = df['g_norm_diff'].max()



# find minimum value of g_norm_ortho - g_norm_meta 
min_g_norm_diff = df['g_norm_diff'].min()

# find average value of g_norm_ortho - g_norm_meta
avg_g_norm_diff = df['g_norm_diff'].mean()

print(F"Maximum E_opt_ortho - E_opt_meta: {max_E_diff} at theta={max_theta}, phi={max_phi}")
print(F"Maximum g_norm_ortho - g_norm_meta: {max_g_norm_diff}")
print(F"Minimum g_norm_ortho - g_norm_meta: {min_g_norm_diff}")
print(F"Average g_norm_ortho - g_norm_meta: {avg_g_norm_diff}")

# print gnorm difference at max_theta and max_phi   
print(F"gnorm difference at max_theta and max_phi: {df[(df['theta'] == max_theta) & (df['phi'] == max_phi)]['g_norm_diff'].values[0]}")

# print ortho gnorm at max_theta and max_phi
print(F"ortho gnorm at max_theta and max_phi: {df[(df['theta'] == max_theta) & (df['phi'] == max_phi)]['g_norm_ortho'].values[0]}")
# print meta gnorm at max_theta and max_phi
print(F"meta gnorm at max_theta and max_phi: {df[(df['theta'] == max_theta) & (df['phi'] == max_phi)]['g_norm_meta'].values[0]}")