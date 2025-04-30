# A Physics-Informed Neural Network for Modeling Higher Order Hydrodynamic Interactions in Heterogeneous Suspensions
# University of British Columbia, August 2024.
# Layal Jbara
#
# This file, Normalize.py, contains functions to normalize the hydrodynamic coefficients in suspension data

#####################################################################################################

"""
Normalize.py contains functions that normalize hydrodynamic coefficients (Cd, Cl, Ct) in suspension data 
based on the volume fraction (phi) of particles. The normalization factor is computed using the cdav function.
"""

#####################################################################################################

import sys
import pandas as pd
import glob
import os

# Add the custom path to the system path for package access
sys.path.append('/zfs/users/lmj44/lmj44/anaconda3/lib/python3.8/site-packages')

# Constant for Cd normalization
Cds = 2.74931

def cdav(phi, a=3.23, b=7.37):
    """Calculate the normalization factor based on the volume fraction (phi)."""
    return a * (1 - phi) ** (-b)

def normalize_cd(df, cd_columns):
    """Normalize Cd columns based on the phi_l column."""
    df['Cd0'] = df['Cd0'].div(cdav(df['phi_l']), axis=0)
    df[cd_columns] = (df[cd_columns] - Cds).div(cdav(df['phi_l']), axis=0)
    return df

def normalize_cl_ct(df, cl_columns, ct_columns):
    """Normalize the Cl and Ct columns based on the phi_l column."""
    df[cl_columns] = df[cl_columns].div(cdav(df['phi_l']), axis=0)
    df[ct_columns] = df[ct_columns].div(cdav(df['phi_l']), axis=0)
    return df

def process_file(file):
    """Process a single file to normalize its data."""
    df = pd.read_csv(file)
    
    # Extract phi value from the file path
    phi_value = file.split('/')[3]
    
    # Identify columns starting with 'Cd', excluding 'Cd0'
    cd_columns = [col for col in df.columns if col.startswith('Cd') and col != 'Cd0']
    df = normalize_cd(df, cd_columns)

    # Identify columns starting with 'Cl' and 'Ct'
    cl_columns = [col for col in df.columns if col.startswith('Cl')]
    ct_columns = [col for col in df.columns if col.startswith('Ct')]

    # Normalize the Cd, Cl, and Ct columns based on the phi_l column
    df = normalize_cl_ct(df, cl_columns, ct_columns)

    # Save the normalized results to a new CSV file
    output_dir = os.path.join('../results/Suspensions/', phi_value)
    os.makedirs(output_dir, exist_ok=True)  # Ensure the directory exists
    df.to_csv(os.path.join(output_dir, 'normalized_results.csv'))
    
def normalize_datasets():
    for file in glob.iglob('../results/Suspensions/*/results.csv'):
        process_file(file)

if __name__ == "__main__":
    normalize_datasets()