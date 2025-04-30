# A Physics-Informed Neural Network for Modeling Higher Order Hydrodynamic Interactions in Heterogeneous Suspensions
# University of British Columbia, August 2024.
# Layal Jbara
#
# This file, Boundary_Neighbour_Search.py, includes preprocessing functions that account for periodic
# boundary conditions when considering the particles.

#####################################################################################################

"""
Boundary_Neighbour_Search.py contains preprocessing functions that search for N nearest neighbors, 
considering periodic boundary conditions.
"""

#####################################################################################################

# Import necessary libraries
import sys
import numpy as np
import polars as pl
import pandas as pd
import re
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import glob
import multiprocessing as mp

# Add custom Python module path
sys.path.append('/zfs/users/lmj44/lmj44/anaconda3/lib/python3.8/site-packages')

def model_function(phi, C=3.470161, p=-7.271069):
    """
    This function calculates the average drag force, \langle F^x \rangle, in a homogeneous suspension
    based on a power law relationship with the volume fraction \langle \phi \rangle.

    The relationship is defined as:
    \langle F^x \rangle = C * (1 - \langle \phi \rangle)^p

    Parameters:
    -----------
    phi: float or array-like.
        The local volume fraction \langle \phi \rangle.
    C: float, optional. 
        The constant. Default is 3.2303.
    p: float, optional. 
        The exponent. Default is -7.3716.

    Returns:
    --------
    float or array.
        The calculated average drag force \langle F^x \rangle.
    """
    return C * (1 - phi) ** p

def find_nearest_neighbors(particles, num_particles, k=50):
    """
    This function computes the N nearest neighbors for each particle in a dataset, 
    represented by a DataFrame. It creates a KDTree for efficient neighbor 
    searching and assigns the nearest neighbor information (indices, coordinates, and distances) 
    to the DataFrame.

    Parameters:
    -----------
    particles: DataFrame.
        A DataFrame containing particle positions, with 'x0' and 'y0' columns.
    num_particles: int. 
        The number of particles to consider.
    k: int, optional. 
        The number of nearest neighbors to find. Default is 50.

    Returns:
    --------
    particles: DataFrame. 
        The original DataFrame with additional columns for nearest neighbors' indices, 
      coordinates, and distances.
    """
    # Create a KDTree for efficient nearest neighbor search
    tree = cKDTree(particles[['x0', 'y0']])

    # Query the nearest neighbors for the first num_particles rows
    distances, indices = tree.query(particles[['x0', 'y0']].iloc[:num_particles], k=k+1)

    # Exclude the point itself by adjusting the indices and distances
    neighbors_indices = indices[:, 1:]
    neighbors_distances = distances[:, 1:]

    # Initialize arrays to store neighbor information
    neighbor_indices = np.zeros((num_particles, k), dtype=int)
    x_values = np.zeros((num_particles, k))
    y_values = np.zeros((num_particles, k))
    r_values = np.zeros((num_particles, k))

    # Populate the arrays with neighbor data
    for i in range(num_particles):
        neighbor_indices[i] = neighbors_indices[i]
        x_values[i] = particles.loc[neighbors_indices[i], 'x0'].values
        y_values[i] = particles.loc[neighbors_indices[i], 'y0'].values
        r_values[i] = neighbors_distances[i]

    # Define column names for the new neighbor information
    neighbor_columns = [f'n{j+1}' for j in range(k)]
    x_columns = [f'x{j+1}' for j in range(k)]
    y_columns = [f'y{j+1}' for j in range(k)]
    r_columns = [f'r{j+1}' for j in range(k)]

    # Assign the neighbor data to the DataFrame
    particles.loc[particles.index[:num_particles], neighbor_columns] = neighbor_indices
    particles.loc[particles.index[:num_particles], x_columns] = x_values
    particles.loc[particles.index[:num_particles], y_columns] = y_values
    particles.loc[particles.index[:num_particles], r_columns] = r_values

    return particles


def mirror(particles, box_size, xmax, ymax, xmin, ymin):
    """
    Mirrors particles near the boundaries and corners of a box, considering periodic boundary conditions.
    
    Parameters:
    ----------
    particles: DataFrame.
        Contains particle coordinates.
    box_size: float, 
        The size of the box used for mirroring.
    xmax, ymax, xmin, ymin: float. 
        The dimensions of the box (max and min in both x and y directions).
    
    Returns:
    --------
    mirrored_particles: DataFrame.
        Contains original and mirrored particle positions.
    """
    # Create a copy of the original particles DataFrame
    mirrored_particles = particles.copy()
    
    # Function to mirror particles and add them to the mirrors list
    def add_mirrored_particles(indices, x_shift=0, y_shift=0):
        mirrored = particles[indices].copy()
        mirrored['x0'] += x_shift
        mirrored['y0'] += y_shift
        return mirrored
    
    mirrors = []
    
    # Mirror particles near boundaries
    if (right_mirror := particles['x0'] > xmax - box_size).any():
        mirrors.append(add_mirrored_particles(right_mirror, x_shift=-(xmax - xmin)))
    
    if (left_mirror := particles['x0'] < xmin + box_size).any():
        mirrors.append(add_mirrored_particles(left_mirror, x_shift=(xmax - xmin)))
    
    if (top_mirror := particles['y0'] > ymax - box_size).any():
        mirrors.append(add_mirrored_particles(top_mirror, y_shift=-(ymax - ymin)))
    
    if (bottom_mirror := particles['y0'] < ymin + box_size).any():
        mirrors.append(add_mirrored_particles(bottom_mirror, y_shift=(ymax - ymin)))
    
    # Mirror particles near corners
    corner_shifts = [
        (right_mirror & top_mirror, -(xmax - xmin), -(ymax - ymin)),
        (left_mirror & top_mirror, (xmax - xmin), -(ymax - ymin)),
        (right_mirror & bottom_mirror, -(xmax - xmin), (ymax - ymin)),
        (left_mirror & bottom_mirror, (xmax - xmin), (ymax - ymin)),
    ]
    
    for mirror_condition, x_shift, y_shift in corner_shifts:
        if mirror_condition.any():
            mirrors.append(add_mirrored_particles(mirror_condition, x_shift, y_shift))
    
    # Concatenate all mirrored DataFrames
    if mirrors:
        mirrored_particles = pd.concat([mirrored_particles] + mirrors)
    
    # Remove duplicates (particles near multiple boundaries or corners)
    mirrored_particles.drop_duplicates(inplace=True)
    
    # Reset index of the concatenated DataFrame
    mirrored_particles.reset_index(drop=True, inplace=True)
    
    return mirrored_particles


def preprocess_data(args):
    """
    This function reads and preprocesses data from a specified test case, including extracting
    diameter, density, and viscosity values from XML and MAC files, computing velocity, and 
    normalizing distances. It computes the hydrodynamic coefficients and concatenates the data
    with porosity information. Additionally, it computes the N nearest neighbors for each particle
    in the dataset and creates a DataFrame containing the reference particle and its neighbors' information.

    Parameters:
    ----------
    args: tuple.
        Contains:
        test_case_name: str.
            The name of the test case directory containing the necessary input files.
        k: int.
            The number of nearest neighbors to consider for each particle.
        Mirroring: bool.
            Whether to apply periodic mirroring.

    Returns:
    --------
    new_df: DataFrame. 
        Preprocessed DataFrame containing computed hydrodynamic coefficients,
        normalized distances, and nearest neighbor information.
    dfp: DataFrame. 
        DataFrame containing the porosity and other related data.
    """
    
    # Unpack arguments
    test_case_name, k, Mirroring = args
    
    # Read diameter value
    radius = extract_value_from_xml(test_case_name + '/Grains/Init/insert.xml', r'<Disc\s+Radius="([^"]+)"\s*/>')
    d = 2. * float(radius)
    
    # Read density and viscosity values
    content = read_file(test_case_name + '/InputFiles/prob_def.mac')
    rho = extract_value_from_pattern(content, r'Density\s*=\s*([0-9.]+)')
    mu = extract_value_from_pattern(content, r'Viscosity\s*=\s*([0-9.]+)')
    
    # Read domain size
    xmin, xmax, ymin, ymax = read_domain_size(content, d)
    
    # Compute velocity
    Re = 10
    U = (Re * mu) / (rho * d)
    
    # Read dataset
    dataset = pd.read_csv(test_case_name + '/DS_results/particle_forces.csv')
    dataset.columns = ['t', 'parID', 'x0', 'y0', 'z0', 'vx0', 'vy0', 'vz0', 'wx0', 'wy0', 'wz0',
                       'Fpx0', 'Fpy0', 'Fpz0', 'Fvx0', 'Fvy0', 'Fvz0', 'Tpx0', 'Tpy0', 'Tpz0',
                       'Tvx0', 'Tvy0', 'Tvz0', 'Tgrad0']
    
    # Extract last timestep data and normalize distances
    df = dataset[dataset['t'] == dataset['t'].max()].reset_index(drop=True)
    df[['x0', 'y0']] /= d
    
    # Compute hydrodynamic coefficients
    df['Cd0'] = (df['Fvx0'] + df['Fpx0']) / (0.5 * rho * U**2 * d)
    df['Cl0'] = (df['Fvy0'] + df['Fpy0']) / (0.5 * rho * U**2 * d)
    df['Ct0'] = (df['Tvz0'] + df['Tpz0']) / (0.5 * rho * U**2 * d**2)
    
    # Concatenate porosity data and compute local volume fraction
    porosity_file = pd.read_csv(test_case_name + '/Res/porosity.csv')
    dfp = pd.concat([df, porosity_file], axis=1)
    dfp['phi_l'] = 1 - dfp['porosity']
    dfp['Cdbar'] = model_function(dfp['phi_l'])
    
    N_p = len(df)
    
    if Mirroring:
        # Compute box side for periodic mirroring
        mini_box_side = compute_mirroring_box(d, k, N_p, xmin, xmax, ymin, ymax, radius)
        new_df = mirror(dfp, mini_box_side, xmax, ymax, xmin, ymin)
    else:
        new_df = dfp
        
    # Find nearest neighbors of unmirrored particles
    new_df = find_nearest_neighbors(new_df, N_p, k)
    
    # Save to CSV
    test_case_name = test_case_name.split('/')[3]
    new_df.head(len(dfp)).to_csv(f'../results/Suspensions/{test_case_name}/boundary_neighbours.csv')

    return new_df, dfp


def read_file(file_path):
    """Reads the contents of a file."""
    with open(file_path, 'r') as file:
        return file.read()


def extract_value_from_xml(file_path, pattern):
    """Extracts a value using a regex pattern from an XML file."""
    content = read_file(file_path)
    match = re.search(pattern, content)
    return match.group(1)


def extract_value_from_pattern(content, pattern):
    """Extracts a value using a regex pattern from the content of a file."""
    match = re.search(pattern, content)
    return float(match.group(1))


def read_domain_size(content, d):
    """Extracts and returns the domain size values, normalized by diameter."""
    patterns = [r'\$DS_Xmin\s*=\s*([0-9.]+)', r'\$DS_Xmax\s*=\s*([0-9.]+)', 
                r'\$DS_Ymin\s*=\s*([0-9.]+)', r'\$DS_Ymax\s*=\s*([0-9.]+)']
    xmin, xmax, ymin, ymax = [float(re.search(pattern, content).group(1)) / d for pattern in patterns]
    return xmin, xmax, ymin, ymax


def compute_mirroring_box(d, k, N_p, xmin, xmax, ymin, ymax, radius):
    """Computes the box side for periodic mirroring."""
    A_p = np.pi * (float(radius) / d)**2
    mini_N_p = 2. * k
    A_domain = (xmax - xmin) * (ymax - ymin)
    A_mini = A_domain * mini_N_p / N_p
    return np.sqrt(A_mini)


def plot_particles(df_new, df_old, d):
    """
    This function plots the positions of particles from two DataFrames, df_new and df_old, with 
    separate colors for the first part and the remaining rows.

    Parameters:
    ----------
    df_new: DataFrame
        Contains the new particle positions.
    df_old: DataFrame 
        Contains the old particle positions (for determining the number of rows to plot in blue).
    d: float.
        Scaling factor for the positions.
    """
    # Create a figure for the plot
    plt.figure(figsize=(12, 12))  # Adjust the width and height as needed

    # Determine the number of rows to plot in blue (based on df_old)
    rows_to_plot_in_blue = len(df_old)

    # Plot the first 'rows_to_plot_in_blue' particles in blue
    plt.scatter(df_new['x0'].iloc[:rows_to_plot_in_blue] * d, 
                df_new['y0'].iloc[:rows_to_plot_in_blue] * d, 
                color='blue', label='First 1000', alpha=0.6)

    # Plot the remaining particles in red
    plt.scatter(df_new['x0'].iloc[rows_to_plot_in_blue:] * d, 
                df_new['y0'].iloc[rows_to_plot_in_blue:] * d, 
                color='red', label='Remaining', alpha=0.6)

    # Adjust plot margins if needed to prevent overlap
    plt.margins(0.1)
    
    # Add labels and legend
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()

    # Show the plot
    plt.show()
    

def plot_neighbours(df, index, k, d):
    """
    This function plots a reference particle and its nearest neighbors from a DataFrame.
    The reference particle is plotted in blue, and its neighbors are plotted in green or red 
    based on a condition.

    Parameters:
    -----------
    df: DataFrame 
        Contains particle data.
    index: int
        The index of the reference particle.
    k: int
        The number of nearest neighbors to plot.
    d: float
        Scaling factor for the positions.
    """
    # Get the reference particle's position
    reference_particle = df.iloc[index]

    # Create a figure for the plot
    plt.figure(figsize=(12, 12))  # Adjust the width and height as needed

    # Plot the reference particle in blue
    plt.scatter(reference_particle['x0'] * d, reference_particle['y0'] * d, color='blue', label='Reference')

    # Plot each of the nearest neighbors
    for i in range(1, k + 1):
        x_col = f'x{i}'
        y_col = f'y{i}'

        # Determine the color based on the condition
        color = 'red' if reference_particle[f'n{i}'] > 9999 else 'green'

        # Plot each neighbor
        plt.scatter(reference_particle[x_col] * d, reference_particle[y_col] * d, color=color)

        # Annotate each neighbor with its index
        plt.annotate(i, 
                     (reference_particle[x_col] * d, reference_particle[y_col] * d), 
                     textcoords="offset points", 
                     xytext=(0, 10), 
                     ha='center', 
                     fontsize=8)

    # Adjust plot margins if needed to prevent overlap
    plt.margins(0.1)

    # Add labels, legend, and show the plot
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.show()
    
    
    

def find_periodic_nearest_neighbours():
    """
    Main function to preprocess all datasets in parallel.
    """
    # Get all dataset paths
    datasets = glob.iglob('../data/Suspensions/*')  

    # Prepare tasks for parallel processing
    tasks = [(file, 45, False) for file in datasets]
    
    # Parallel mode preprocessing using multiprocessing
    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.map(preprocess_data, tasks)
        
        
if __name__ == "__main__":
    find_periodic_nearest_neighbours()