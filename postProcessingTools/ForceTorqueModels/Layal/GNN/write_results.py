"""
Graph Neural Networks for Modeling Hydrodynamic Closure Laws in Non-Spherical Particle-Laden Flows
University of British Columbia, March 2025
Author: Layal Jbara

This script, write_results.py:
- Extracts the maximum validation R² scores from GNN training results across multiple simulation cases 
  and compiles them into a summary CSV file for further analysis.

"""

import os
import json
import pandas as pd

# -------------------------------
# Configuration Parameters
# -------------------------------
shapes = ['Tetra', 'Cube', 'Octa', 'Dodeca', 'Icosa']
Phis = ['Phi005', 'Phi01', 'Phi02']
Res = ['Re1', 'Re10', 'Re100']
neighbours = ['10', '20', '30', '40', '50']

# Initialize an empty list to store results
results = []

# -------------------------------
# Loop through all combinations
# -------------------------------
for phi in Phis:
    for Re in Res:
        for shape in shapes:
            for N in neighbours:
                # Construct the path to the JSON result file
                json_file = f'../results/{phi}/{Re}/{shape}/Fx_det_results_84channels_neighbours_{N}_BS_size_128.json'

                # Check if the file exists before trying to open
                if os.path.exists(json_file):
                    try:
                        # Load JSON file contents
                        with open(json_file, 'r') as file:
                            data = json.load(file)

                        # Validate that data is a list of dictionaries
                        if isinstance(data, list):
                            # Extract the maximum "best_r2_valid" value in the file
                            max_r2_valid = max(
                                (entry.get("best_r2_valid", float('-inf')) for entry in data),
                                default=float('-inf')
                            )

                            # Only add entry if a valid max was found
                            if max_r2_valid > float('-inf'):
                                results.append({
                                    "Phi": phi,
                                    "Re": Re,
                                    "Shape": shape,
                                    "Neighbours": N,
                                    "Max_best_r2_valid": max_r2_valid
                                })
                    except json.JSONDecodeError:
                        print(f"[ERROR] Failed to decode JSON: {json_file}")
                else:
                    print(f"[WARNING] File not found: {json_file}")

# -------------------------------
# Convert Results to DataFrame
# -------------------------------
df = pd.DataFrame(results)

# Save to CSV file
output_csv = "../results/Fx_best_r2_results.csv"
df.to_csv(output_csv, index=False)
print(f"[INFO] Results saved to {output_csv}")



