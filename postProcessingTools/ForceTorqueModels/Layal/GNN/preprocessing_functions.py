# Import necessary libraries
#######################################################################
import os
import numpy as np
import pandas as pd
import re
from scipy.spatial import cKDTree
import glob
import multiprocessing as mp
from scipy.interpolate import RegularGridInterpolator
#######################################################################


def mirror_3D(particles, box_size, xmax, ymax, xmin, ymin, zmax, zmin):
    """
    Mirrors particles near the boundaries and corners of a box, considering periodic boundary conditions.
    
    Parameters:
    ----------
    particles: DataFrame
        Contains particle coordinates.
    box_size: float
        The size of the box used for mirroring.
    xmax, ymax, zmax, xmin, ymin, zmin: float
        Maximum and minimum dimensions of the box in the x, y, and z directions.
    
    Returns:
    --------
    mirrored_particles: DataFrame
        Contains original and mirrored particle positions.
    """
    
     # Create a copy of the original particles DataFrame
    mirrored_particles = particles.copy()
    
    # Calculate mirrored positions for each boundary and corner
    mirrors = []

    # Function to mirror particles
    def add_mirrored_particles(indices, x_shift=0, y_shift=0, z_shift=0):
        mirrored = particles[indices].copy()
        mirrored['x0'] += x_shift
        mirrored['y0'] += y_shift
        mirrored['z0'] += z_shift
        return mirrored

    # Mirror particles near faces (1D shifts)
    if (right_mirror := particles['x0'] > xmax - box_size).any():
        mirrors.append(add_mirrored_particles(right_mirror, x_shift=-(xmax - xmin)))

    if (left_mirror := particles['x0'] < xmin + box_size).any():
        mirrors.append(add_mirrored_particles(left_mirror, x_shift=(xmax - xmin)))

    if (top_mirror := particles['y0'] > ymax - box_size).any():
        mirrors.append(add_mirrored_particles(top_mirror, y_shift=-(ymax - ymin)))

    if (bottom_mirror := particles['y0'] < ymin + box_size).any():
        mirrors.append(add_mirrored_particles(bottom_mirror, y_shift=(ymax - ymin)))

    if (front_mirror := particles['z0'] > zmax - box_size).any():
        mirrors.append(add_mirrored_particles(front_mirror, z_shift=-(zmax - zmin)))

    if (back_mirror := particles['z0'] < zmin + box_size).any():
        mirrors.append(add_mirrored_particles(back_mirror, z_shift=(zmax - zmin)))

    # Mirror particles near edges (2D shifts)
    edge_shifts = [
        (right_mirror & top_mirror, -(xmax - xmin), -(ymax - ymin), 0),
        (left_mirror & top_mirror, (xmax - xmin), -(ymax - ymin), 0),
        (right_mirror & bottom_mirror, -(xmax - xmin), (ymax - ymin), 0),
        (left_mirror & bottom_mirror, (xmax - xmin), (ymax - ymin), 0),
       
        (right_mirror & front_mirror, -(xmax - xmin), 0, -(zmax - zmin)),
        (left_mirror & front_mirror, (xmax - xmin), 0, -(zmax - zmin)),
        (right_mirror & back_mirror, -(xmax - xmin), 0, (zmax - zmin)),
        (left_mirror & back_mirror, (xmax - xmin), 0, (zmax - zmin)),
        
#         (top_mirror & front_mirror, 0, -(ymax - ymin), -(zmax - zmin)),
#         (bottom_mirror & front_mirror, 0, (ymax - ymin), -(zmax - zmin)),
#         (top_mirror & back_mirror, 0, -(ymax - ymin), (zmax - zmin)),
#         (bottom_mirror & back_mirror, 0, (ymax - ymin), (zmax - zmin)),
    ]

    # Add mirrored particles for edges
    for mirror_condition, x_shift, y_shift, z_shift in edge_shifts:
        if mirror_condition.any():
            mirrors.append(add_mirrored_particles(mirror_condition, x_shift, y_shift, z_shift))

#     # Mirror particles near corners (3D shifts)
#     corner_shifts = [
#        (right_mirror & top_mirror & front_mirror, -(xmax - xmin), -(ymax - ymin), -(zmax - zmin)),
#        (left_mirror & top_mirror & front_mirror, (xmax - xmin), -(ymax - ymin), -(zmax - zmin)),
#        (right_mirror & bottom_mirror & front_mirror, -(xmax - xmin), (ymax - ymin), -(zmax - zmin)),
#        (left_mirror & bottom_mirror & front_mirror, (xmax - xmin), (ymax - ymin), -(zmax - zmin)),
#        (right_mirror & top_mirror & back_mirror, -(xmax - xmin), -(ymax - ymin), (zmax - zmin)),
#        (left_mirror & top_mirror & back_mirror, (xmax - xmin), -(ymax - ymin), (zmax - zmin)),
#        (right_mirror & bottom_mirror & back_mirror, -(xmax - xmin), (ymax - ymin), (zmax - zmin)),
#        (left_mirror & bottom_mirror & back_mirror, (xmax - xmin), (ymax - ymin), (zmax - zmin)),
#    ]

#     # Add mirrored particles for corners
#     for mirror_condition, x_shift, y_shift, z_shift in corner_shifts:
#         if mirror_condition.any():
#             mirrors.append(add_mirrored_particles(mirror_condition, x_shift, y_shift, z_shift))

    # Concatenate all mirrored DataFrames into a single DataFrame
    if mirrors:
        mirrored_particles = pd.concat([mirrored_particles] + mirrors)
    
    # Remove duplicate particles
    mirrored_particles.drop_duplicates(inplace=True)
    
    # Reset the index of the resulting DataFrame
    mirrored_particles.reset_index(drop=True, inplace=True)
    
    return mirrored_particles



def find_nearest_neighbors(particles, num_particles, k=50):
    """
    Finds the nearest neighbors for a given set of particles and calculates 
    relevant distances and spherical angles. This function uses a KDTree 
    structure for efficient nearest-neighbor search.

    Parameters:
    ----------
    particles : DataFrame
        The DataFrame containing the particle data, which must include 'parID', 'x0', 'y0', and 'z0' columns 
        representing the particle ID and their 3D coordinates.
    num_particles : int
        The number of particles to consider for neighbor search.
    k : int, optional, default=50
        The number of nearest neighbors to find for each particle. The default is 50.

    Returns:
    --------
    particles : DataFrame
        The original DataFrame with additional columns for the nearest neighbors' 
        data (coordinates, relative distances, and spherical angles).
    """
    
    # Create a KDTree with 3D coordinates for efficient nearest-neighbor search
    tree = cKDTree(particles[['x0', 'y0', 'z0']])
    
    # Query the nearest neighbors for the first `num_particles` particles
    distances, indices = tree.query(particles[['x0', 'y0', 'z0']].iloc[:num_particles], k=k+1)
    
    # Exclude the particle itself from the neighbors
    neighbors_indices = indices[:, 1:]  # Exclude the first column (self)
    neighbors_distances = distances[:, 1:]  # Exclude the first column (self)

    # Initialize arrays for storing neighbors' data
    neighbor_keys = ['n', 'x', 'y', 'z', 'r', 'rx', 'ry', 'rz', 'theta', 'phi']
    neighbor_data = {key: np.zeros((num_particles, k), dtype=int if key == 'n' else float) 
                     for key in neighbor_keys}    

    # Loop over each particle and assign data for its neighbors
    for i in range(num_particles):
        # Get the parID and coordinates of neighbors for the current particle
        neighbor_data['n'][i] = particles.loc[neighbors_indices[i], 'parID'].values
        neighbor_data['x'][i] = particles.loc[neighbors_indices[i], 'x0'].values
        neighbor_data['y'][i] = particles.loc[neighbors_indices[i], 'y0'].values
        neighbor_data['z'][i] = particles.loc[neighbors_indices[i], 'z0'].values
        neighbor_data['r'][i] = neighbors_distances[i]
        
        # Calculate relative distances from the current particle to each of its neighbors
        particle_index = particles.index[i]
        neighbor_data['rx'][i] = neighbor_data['x'][i] - particles.loc[particle_index, 'x0']
        neighbor_data['ry'][i] = neighbor_data['y'][i] - particles.loc[particle_index, 'y0']
        neighbor_data['rz'][i] = neighbor_data['z'][i] - particles.loc[particle_index, 'z0']
        
        # Calculate spherical angles (theta and phi) for each neighbor
        neighbor_data['theta'][i] = np.arccos(neighbor_data['rz'][i] / neighbor_data['r'][i])
        neighbor_data['phi'][i] = np.arctan2(neighbor_data['ry'][i], neighbor_data['rx'][i])

    # Assign the calculated neighbor data arrays to the DataFrame
    for key, values in neighbor_data.items():
        # Create column names for each neighbor (e.g., n1, n2, ..., n50)
        cols = [f'{key}{j+1}' for j in range(k)]
        particles.loc[particles.index[:num_particles], cols] = values

    # Return the DataFrame with added neighbor data
    return particles


def reflect_data(df, flip_y=False, flip_z=False):
    """
    This function prepares a copy of the input DataFrame by optionally reflecting the values of 
    the 'y' and 'z' columns based on the specified flags. The flipping operation multiplies 
    the values of the relevant columns by -1.

    Parameters:
    - df: DataFrame
        The input DataFrame containing the data to be processed.
    - flip_y: bool, optional, default=False
        If True, flips all columns starting with 'y'.
    - flip_z: bool, optional, default=False
        If True, flips all columns starting with 'z'.

    Returns:
    - new_df: DataFrame
        A new DataFrame with flipped values in the 'y' and/or 'z' columns, if specified.
    """
    new_df = df.copy()

    if flip_y:
        for col in new_df.columns:
            if col.startswith('y'):
                new_df[col] = -df[col]

    if flip_z:
        for col in new_df.columns:
            if col.startswith('z'):
                new_df[col] = -df[col]

    return new_df


def get_two_angles(M, resolution=100):
    import numpy as np
    
    phi = np.linspace(0, 2. * np.pi, resolution)
    psi = np.linspace(0, 2. * np.pi, resolution)
    
    _phi, _psi = np.meshgrid(phi, psi)
    _residual = np.zeros((resolution, resolution))

    span_dir_1 = np.array([0, 1, 0])
    span_dir_2 = np.array([0, 0, 1])
    
    for i in range(resolution):
        for j in range(resolution):
            R1 = Rodrigues_Rotation(span_dir_2, _phi[i, j])
            p = np.matmul(R1, span_dir_1)
            R2 = Rodrigues_Rotation(p, _psi[i, j])
            _M = np.matmul(R2, R1)
            diff = np.subtract(_M, M)
            _residual[i, j] = np.linalg.norm(diff)
    
    min_i = np.argmin(_residual)
    return _phi.flatten()[min_i], _psi.flatten()[min_i]


def Rodrigues_Rotation(vec, theta):
    import numpy as np
    
    R = np.ones((3, 3))
    
    R[0, 0] = np.cos(theta) + vec[0]**2 * (1. - np.cos(theta))
    R[0, 1] = vec[0] * vec[1] * (1. - np.cos(theta)) - vec[2] * np.sin(theta)
    R[0, 2] = vec[1] * np.sin(theta) + vec[0] * vec[2] * (1. - np.cos(theta))
    
    R[1, 0] = vec[2] * np.sin(theta) + vec[0] * vec[1] * (1. - np.cos(theta))
    R[1, 1] = np.cos(theta) + vec[1]**2 * (1. - np.cos(theta))
    R[1, 2] = vec[1] * vec[2] * (1. - np.cos(theta)) - vec[0] * np.sin(theta)
    
    R[2, 0] = vec[0] * vec[2] * (1. - np.cos(theta)) - vec[1] * np.sin(theta)
    R[2, 1] = vec[0] * np.sin(theta) + vec[1] * vec[2] * (1. - np.cos(theta))
    R[2, 2] = np.cos(theta) + vec[2]**2 * (1. - np.cos(theta))
    
    return R


def insertA_reader(file):

    
    RotationM = pd.DataFrame(columns=['parID', 'x', 'y', 'z', 'm11', 'm12', 'm13', 'm21', 'm22', 'm23', 'm31', 'm32', 'm33'])
    
    with open(file, 'r') as f:
        b_read = False
        Nread = 0
        i = 0
        
        for line in f:
            if i < Nread:
                arr = line.split()
                new_row = pd.DataFrame({
                    'parID': int(arr[0]),
                    'x': float(arr[3]),
                    'y': float(arr[4]),
                    'z': float(arr[5]),
                    'm11': float(arr[6]),
                    'm12': float(arr[7]),
                    'm13': float(arr[8]),
                    'm21': float(arr[9]),
                    'm22': float(arr[10]),
                    'm23': float(arr[11]),
                    'm31': float(arr[12]),
                    'm32': float(arr[13]),
                    'm33': float(arr[14]),
                }, index=[0])
                RotationM = pd.concat([new_row, RotationM.loc[:]]).reset_index(drop=True)
                i += 1
            
            if 'Texte' in line:
                Nread = int(line.split()[1])
                b_read = True
    
    RotationM = RotationM.sort_values(by=['parID']).reset_index(drop=True)
    return RotationM




def process_angle(test_case_name,shape,reflect_y=False,reflect_z=False):
    
    if shape == "Cube" or shape=='Cube2' or shape =='Cube3' or shape=='Cube4':
        refshape = "cube"
    elif shape == "Tetra":
        refshape = "isotetra"
    elif shape == "Octa":
        refshape = "octaedre"
    elif shape == "Dodeca":
        refshape = "dodecaedre"
    elif shape == "Icosa":
        refshape = "icosaedre"
    
    
    rotation_file_path = os.path.join(test_case_name, 'Grains', 'Init', 'insertA.result')
    RotationM = insertA_reader(rotation_file_path)
    ref_file = os.path.join(test_case_name, 'Grains', 'Init', f'{refshape}.insert')

    
    # Reading the number of vertices and location for reference Platonic
    Nvertex = 0
    vertex = []
    with open(ref_file,'r') as f:
        i = 0
        for line in f:
            if (i == 1):
                Nvertex = int(line.split()[0])
            if ((i > 1) and (i <= Nvertex+1)):
                vertex.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])] )
            i = i + 1

    def Rotate(vector,M):
        vx = vector[0]*M[0,0] + vector[1]*M[0,1] + vector[2]*M[0,2];
        vy = vector[0]*M[1,0] + vector[1]*M[1,1] + vector[2]*M[1,2];
        vz = vector[0]*M[2,0] + vector[1]*M[2,1] + vector[2]*M[2,2];

        return np.array([vx, vy, vz])

    VertexDataFrame = pd.DataFrame(columns=['parID', 'x', 'y', 'z'])

    for i in range(len(RotationM)):
        MRot = np.array([[RotationM.m11.iloc[i], RotationM.m12.iloc[i], RotationM.m13.iloc[i]],
                         [RotationM.m21.iloc[i], RotationM.m22.iloc[i], RotationM.m23.iloc[i]],
                         [RotationM.m31.iloc[i], RotationM.m32.iloc[i], RotationM.m33.iloc[i]]])

        for j in range(Nvertex):
            v_org = np.zeros(3)
            v_org[0] = vertex[j][0]
            v_org[1] = vertex[j][1]
            v_org[2] = vertex[j][2]

            v_trans = Rotate(v_org,MRot)

            v_trans[0] = v_trans[0] + RotationM.x.iloc[i]
            v_trans[1] = v_trans[1] + RotationM.y.iloc[i]
            v_trans[2] = v_trans[2] + RotationM.z.iloc[i]

            new_row = pd.DataFrame({'parID':RotationM.parID.iloc[i], 
                                    'x':v_trans[0],
                                    'y':v_trans[1],
                                    'z':v_trans[2],},
                                   index=[0])
            VertexDataFrame = pd.concat([new_row,VertexDataFrame.loc[:]]).reset_index(drop=True)

    VertexDataFrame = VertexDataFrame.sort_values(by=['parID']).reset_index(drop=True)


    two_angles = pd.DataFrame(columns=['parID', 'x', 'y', 'z', 'phi', 'psi'])

    for i in range(len(RotationM)):
        M_GRAINS = np.array([[RotationM.m11.iloc[i], RotationM.m12.iloc[i], RotationM.m13.iloc[i]],
                             [RotationM.m21.iloc[i], RotationM.m22.iloc[i], RotationM.m23.iloc[i]],
                             [RotationM.m31.iloc[i], RotationM.m32.iloc[i], RotationM.m33.iloc[i]]])

        print(f"Calculating phi and psi for {i}")
        
        if reflect_y==True and reflect_z==False:
        # Define the reflection matrix for flipping in the y-direction
            R_reflect_y = np.array([
                [1,  0,  0],
                [0, -1,  0],
                [0,  0,  1]
            ])

            # Apply the reflection transformation
            M_GRAINS_reflected = R_reflect_y @ M_GRAINS @ R_reflect_y

        
        elif reflect_y==False and reflect_z==True:
            R_reflect_z = np.array([
                    [1,  0,  0],
                    [0,  1,  0],
                    [0,  0, -1]
                ])
            # Apply the reflection transformation
            M_GRAINS_reflected = R_reflect_z @ M_GRAINS @ R_reflect_z

            
        elif reflect_y==True and reflect_z==True:
            # Define the reflection matrix for flipping in the y-direction
            R_reflect_yz = np.array([
                [1,  0,  0],
                [0, -1,  0],
                [0,  0, -1]
            ])


            # Apply the reflection transformation
            M_GRAINS_reflected = R_reflect_yz @ M_GRAINS @ R_reflect_yz

        else:
            M_GRAINS_reflected=M_GRAINS
            
        phi, psi = get_two_angles(M_GRAINS_reflected)
        new_row = pd.DataFrame({'parID':RotationM.parID.iloc[i], 
                                'x':RotationM.x.iloc[i],
                                'y':RotationM.y.iloc[i],
                                'z':RotationM.z.iloc[i],
                                'phi':phi % (np.pi/2.),
                                'psi':psi % (np.pi/2.),}, 
                                 index=[0])
        two_angles = pd.concat([new_row,two_angles.loc[:]]).reset_index(drop=True)

    two_angles=two_angles.sort_values(by=['parID'])
    two_angles = two_angles.reset_index(drop=True)


    if shape == 'Octa' and phi=='Phi005':
        df_dropped = two_angles[two_angles['parID'] != 706].copy()  
        df_dropped.loc[:, 'parID'] = range(0, len(df_dropped) )
        df_dropped = df_dropped.reset_index(drop=True)
        two_angles=df_dropped 
   
    return    two_angles 
