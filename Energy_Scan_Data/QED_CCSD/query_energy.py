import pandas as pd
import numpy as np
import sys

def load_data(file_name):
    """Loads the data and cleans it to handle headers and string rows."""
    try:
        # Load, skipping row 1 (the dashed line)
        df = pd.read_csv(file_name, delim_whitespace=True, skiprows=[1])
        
        # Force numeric, turning non-numeric junk into NaN, then drop NaNs
        numeric_cols = ['theta', 'phi', 'Para_E', 'Ortho_E', 'Meta_E']
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        df = df.dropna(subset=numeric_cols)
        return df
    except FileNotFoundError:
        print(f"Error: {file_name} not found.")
        sys.exit(1)

def find_most_stabilizing_orientation(df):
    """Scans the entire grid to find the most favorable stabilization."""
    AU_TO_KCAL = 627.509
    pairs = [
        ('Ortho_E', 'Meta_E', 'Ortho - Meta'),
        ('Ortho_E', 'Para_E',  'Ortho - Para'),
        ('Para_E',  'Meta_E',  'Para - Meta')
    ]
    
    print(f"\n{'='*50}\nGLOBAL STABILIZATION SEARCH\n{'='*50}")
    
    for col1, col2, name in pairs:
        diff_series = (df[col1] - df[col2]) * AU_TO_KCAL
        min_idx = diff_series.idxmin()
        min_val = diff_series.loc[min_idx]
        row = df.loc[min_idx]
        
        print(f"\nMost favorable {name}:")
        print(f"  Max Stabilization: {min_val:8.4f} kcal/mol")
        print(f"  At Angle:          Theta = {row['theta']:.1f}°, Phi = {row['phi']:.1f}°")

def query_energies(df):
    """Interactive loop to query specific orientations."""
    AU_TO_KCAL = 627.509
    
    print(f"\n{'='*50}\nINTERACTIVE ENERGY QUERY TOOL\n{'='*50}")
    print(f"File loaded with {len(df)} data points.")
    
    while True:
        user_in = input("\nEnter target (theta, phi) in degrees (or 'q' to quit): ")
        if user_in.lower() == 'q': break
        
        try:
            t_target, p_target = map(float, user_in.split(','))
            
            # Nearest neighbor search
            distances = np.sqrt((df['theta'] - t_target)**2 + (df['phi'] - p_target)**2)
            idx = distances.idxmin()
            row = df.loc[idx]
            
            e_p, e_o, e_m = row['Para_E'], row['Ortho_E'], row['Meta_E']
            
            print(f"\nResults for closest point ({row['theta']:.1f}°, {row['phi']:.1f}°):")
            print(f"{'Species':<10} {'Energy (Ha)':<20}")
            print(f"{'-'*30}")
            print(f"{'Para':<10} {e_p:.8f}")
            print(f"{'Ortho':<10} {e_o:.8f}")
            print(f"{'Meta':<10} {e_m:.8f}")
            
            print(f"\nEnergy Differences (kcal/mol):")
            print(f"Ortho - Meta:  {(e_o - e_m) * AU_TO_KCAL:8.4f}")
            print(f"Ortho - Para:  {(e_o - e_p) * AU_TO_KCAL:8.4f}")
            print(f"Para  - Meta:  {(e_p - e_m) * AU_TO_KCAL:8.4f}")
            
        except (ValueError, TypeError):
            print("Invalid input! Please enter two numbers separated by a comma (e.g., 90, 180).")

if __name__ == "__main__":
    FILE_NAME = "CCSD_Combined_Results.txt"
    data = load_data(FILE_NAME)
    
    # 1. Run the global scan summary
    find_most_stabilizing_orientation(data)
    
    # 2. Start the interactive tool
    query_energies(data)
