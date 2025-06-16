# Graph Neural Networks for Modeling Hydrodynamic Closure Laws in Non-Spherical Particle-Laden Flows
# University of British Columbia, March 2025.
# Layal Jbara
#
# This file, Boundary_Neighbour_Search.py, includes preprocessing functions that account for periodic
# boundary conditions when considering the particles.

#####################################################################################################

"""
Preprocess_Data.py

This script contains functions for preprocessing data in the context of particle simulations. The main tasks include:
- Searching for and applying periodic boundary conditions with mirrored data.
- Computing the N nearest neighbors for each particle in the dataset.
- Normalizing the spatial coordinates.
- Calculating hydrodynamic coefficients based on the particle forces and velocities.
- Saving the processed data into a CSV file for further analysis.

"""

#####################################################################################################

# Import necessary libraries
#######################################################################
import os
import numpy as np
import pandas as pd
import re
from scipy.spatial import cKDTree
import glob
import multiprocessing as mp
import preprocessing_functions as prf
#######################################################################


def preprocess_data(args):
    """
    This function reads and preprocesses data from a specified test case, including extracting
    diameter, density, viscosity values from XML and MAC files, computing velocity, and 
    normalizing distances. It also computes the hydrodynamic coefficients, calculates the N 
    nearest neighbors for each particle in the dataset, and creates a DataFrame containing 
    the reference particle and its neighbors' information.

    Args:
    - test_case_name (str): The name of the test case directory containing the necessary input files.
    - k (int): The number of nearest neighbors to consider for each particle.
    - mirror_func (bool): A flag indicating whether or not to apply periodic mirroring to the data.

    Returns:
    - pd.DataFrame: Preprocessed DataFrame containing computed hydrodynamic coefficients, 
                     normalized distances, and nearest neighbor information.
    """
    
    # Unpack the arguments
    test_case_name, k, mirror_func = args
    print(test_case_name)

    # Extract Phi, Re, and Shape from the test_case_name path
    parts = test_case_name.split('/')
    phi = parts[-3]  # e.g., "Phi005"
    Re = parts[-2]   # e.g., "Re1"
    shape = parts[-1]  # e.g., "Cube" or other shape

    # Print extracted values for debugging
    print(phi)
    print(Re)
    print(shape)

    # Create the output directory if it doesn't exist
    output_dir = os.path.join('../results', phi, Re, shape)
    os.makedirs(output_dir, exist_ok=True)
    
    ## Extract Diameter Value (hardcoded for now)
    d = 0.001
    
    # Extract Density and Viscosity Values from the .mac file
    with open(test_case_name + '/InputFiles/prob_def.mac', 'r') as file:
        content = file.read()

    # Extract density value
    pattern_density = r'DS_rho\s*=\s*([0-9.]+)'
    match = re.search(pattern_density, content)
    rho = float(match.group(1))
    
    # Extract viscosity value
    pattern_viscosity = r'DS_mu\s*=\s*([0-9.]+)'
    match = re.search(pattern_viscosity, content)
    mu = float(match.group(1))
    
    # Extract velocity value (Umean)
    pattern_velocity = r'DS_Umean\s*=\s*([0-9.]+)'   
    match = re.search(pattern_velocity, content)
    U = float(match.group(1))
    
    # Extract domain size (xmin, xmax, ymin, ymax, zmin, zmax)
    patterns = {
        'xmin': r'\$DS_Xmin\s*=\s*([0-9.]+)',
        'xmax': r'\$DS_Xmax\s*=\s*([0-9.]+)',
        'ymin': r'\$DS_Ymin\s*=\s*([0-9.]+)',
        'ymax': r'\$DS_Ymax\s*=\s*([0-9.]+)',
        'zmin': r'\$DS_Zmin\s*=\s*([0-9.]+)',
        'zmax': r'\$DS_Zmax\s*=\s*([0-9.]+)',
    }
    domain_size = {key: float(re.search(pattern, content).group(1)) for key, pattern in patterns.items()}
    
    # Normalize domain size by diameter
    xmin, xmax, ymin, ymax, zmin, zmax = (domain_size[key] / d for key in ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'])
    
    # Read the particle force data into a DataFrame
    dataset = pd.read_csv(test_case_name + '/DS_results/particle_forces.csv')
    dataset.columns = ['t', 'parID', 'x0', 'y0', 'z0', 'vx0', 'vy0', 'vz0', 'wx0', 'wy0', 'wz0',
                       'Fpx0', 'Fpy0', 'Fpz0', 'Fvx0', 'Fvy0', 'Fvz0', 'Tpx0', 'Tpy0', 'Tpz0',
                       'Tvx0', 'Tvy0', 'Tvz0', 'Tgrad0']
    
    # Extract the last timestep data
    df = dataset[dataset['t'] == dataset['t'].max()].reset_index(drop=True)

    # Remove specific particle (for a certain shape and phi)
    if shape == 'Octa' and phi == 'Phi005':
        df_dropped = df[df['parID'] != 706].copy()
        df_dropped.loc[:, 'parID'] = range(0, len(df_dropped))
        df_dropped = df_dropped.reset_index(drop=True)
        df = df_dropped
    
    # Normalize particle coordinates by diameter
    df[['x0', 'y0', 'z0']] /= d
    
    F_scale =  3.*np.pi*mu*d*U
    T_scale =  3.*np.pi*mu*d*d*U
    
    ## Compute hydrodynamic coefficients
    df['Cd0']  = (df['Fvx0']+df['Fpx0'])/(F_scale)
    df['Cly0'] = (df['Fvy0']+df['Fpy0'])/(F_scale)
    df['Clz0'] = (df['Fvz0']+df['Fpz0'])/(F_scale)
    df['Ctz0'] = (df['Tvz0']+df['Tpz0'])/(T_scale)
    df['Cty0'] = (df['Tvy0']+df['Tpy0'])/(T_scale)
    df['Ctx0'] = (df['Tvx0']+df['Tpx0'])/(T_scale)
    
    
    # Compute hydrodynamic coefficients (Cd, Cly, Clz, etc.)
#    df['Cd0'] = (df['Fvx0'] + df['Fpx0']) / (0.5 * rho * U**2 * np.pi * (d**2 / 4))
#    df['Cly0'] = (df['Fvy0'] + df['Fpy0']) / (0.5 * rho * U**2 * np.pi * (d**2 / 4))
#    df['Clz0'] = (df['Fvz0'] + df['Fpz0']) / (0.5 * rho * U**2 * np.pi * (d**2 / 4))
#    df['Ctz0'] = (df['Tvz0'] + df['Tpz0']) / (0.5 * rho * U**2 * np.pi * (d**3 / 8))
#    df['Cty0'] = (df['Tvy0'] + df['Tpy0']) / (0.5 * rho * U**2 * np.pi * (d**3 / 8))
#    df['Ctx0'] = (df['Tvx0'] + df['Tpx0']) / (0.5 * rho * U**2 * np.pi * (d**3 / 8))
    
 ##########################################################################################
    # Apply periodic mirroring if specified
    if mirror_func:
        mini_N_p = k
        V_domain = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
        V_mini = V_domain * mini_N_p / len(df)
        mini_box_side = V_mini ** (1/3)
        
        # Apply mirroring
        new_df = prf.mirror_3D(df, mini_box_side, xmax, ymax, xmin, ymin, zmax, zmin)
    else:
        new_df = df
    ##########################################################################################
    # Find nearest neighbors of particles (mirrored or not)
    new_df = prf.find_nearest_neighbors(new_df, len(df), k)
    # Save the processed data to a CSV file
    output_file = os.path.join(output_dir, 'boundary_neighbours_optimized.csv')
    new_df.head(len(df)).to_csv(output_file, index=False)
    ##########################################################################################
     ## Reflection Symmetry    
    data_y_flipped = prf.reflect_data(new_df,  flip_y=True, flip_z=False)
    data_z_flipped = prf.reflect_data(new_df,  flip_y=False, flip_z=True)
    data_zy_flipped = prf.reflect_data(new_df,  flip_y=True, flip_z=True)

    output_file = os.path.join(output_dir, 'boundary_neighbours_flip_y_optimized.csv')
    data_y_flipped.to_csv(output_file, index=False)

    output_file = os.path.join(output_dir, 'boundary_neighbours_flip_z_optimized.csv')
    data_z_flipped.to_csv(output_file, index=False)

    output_file = os.path.join(output_dir, 'boundary_neighbours_flip_zy_optimized.csv')
    data_zy_flipped.to_csv(output_file, index=False)
    ##########################################################################################
    ## Rotation Angles Computation Symmetry    
    angle_data=prf.process_angle(test_case_name,shape,reflect_y=False,reflect_z=False)
    output_file = os.path.join(output_dir, 'angle_optimized.csv')        
    angle_data.to_csv(output_file, index=False)

    angle_data_y=prf.process_angle(test_case_name,shape,reflect_y=True,reflect_z=False)
    output_file = os.path.join(output_dir, 'angle_flip_y_optimized.csv')        
    angle_data_y.to_csv(output_file, index=False)

    angle_data_z=prf.process_angle(test_case_name,shape,reflect_y=False,reflect_z=True)
    output_file = os.path.join(output_dir, 'angle_flip_z_optimized.csv')        
    angle_data_z.to_csv(output_file, index=False)

    angle_data_yz=prf.process_angle(test_case_name,shape,reflect_y=True,reflect_z=True)
    output_file = os.path.join(output_dir, 'angle_flip_yz_optimized.csv')        
    angle_data_yz.to_csv(output_file, index=False)


    
    return new_df

def find_periodic_nearest_neighbours():
    """
    Main function to preprocess all datasets in parallel.
    """
    # Get all dataset paths
    datasets = glob.iglob('../data/*/*/*')  


    # Prepare tasks for parallel processing
    tasks = [(file, 150, True) for file in datasets]
    
    num_cpus = int(os.getenv('SLURM_CPUS_ON_NODE', mp.cpu_count()))
    print(f"Using {num_cpus} CPUs for processing")
    print(tasks)
    # Parallel mode preprocessing using multiprocessing
    with mp.Pool(processes=num_cpus) as pool:
        pool.map(preprocess_data, tasks)



if __name__ == "__main__":
     find_periodic_nearest_neighbours()
