# A Physics-Informed Neural Network for Modeling Higher Order Hydrodynamic Interactions in Heterogeneous Suspensions
# University of British Columbia, August 2024.
# Layal Jbara
#
# This file, Preprocess_Torque.py, contains functions that preprocess raw data, generate features for higher-order interactions, 
# apply symmetry-based augmentation, and prepare labels for training a Torque Hierarchical Physics-Informed Neural Network (PINN). 
# It also handles shuffling and concatenation of data to improve model generalization.
#####################################################################################################

"""
Preprocess_Torque.py contains functions for preprocessing data, generating features, applying symmetry augmentation,
and preparing labels for training a Torque Hierarchical PINN.
"""

#####################################################################################################

import sys
import numpy as np
import pandas as pd
from itertools import combinations
import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
from tensorflow.keras.layers import Concatenate, Dense, Dropout, Input, Lambda
from tensorflow.keras.models import Model
import warnings
import argparse
import random
from sklearn.model_selection import KFold, train_test_split



# Add custom path for module imports
sys.path.append('/zfs/users/lmj44/lmj44/anaconda3/lib/python3.8/site-packages')

# Set logging level and filter warnings
tf.get_logger().setLevel('ERROR')
warnings.filterwarnings('ignore')

# Print the number of available GPUs
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

try:
    # Disable all GPUs
    tf.config.set_visible_devices([], 'GPU')
    visible_devices = tf.config.get_visible_devices()
    
    # Assert that no GPU devices are visible
    for device in visible_devices:
        assert device.device_type != 'GPU'
except Exception as e:
    # Handle exceptions (e.g., invalid device or cannot modify virtual devices)
    print(f"Error occurred: {e}")

#####################################################################################################

def dataset_process(csv_file, N_2, N_3, N_4):
    """
    Processes the raw dataset and generates features for higher-order interactions.

    This function reads a dataset from a CSV file, extracts and transforms the data for higher-order interactions 
    (binary, ternary, and quaternary). The extracted features are reshaped into arrays suitable for model training.
    It also calculates relative distances for each partiCte pair and higher-order combinations, and returns these 
    transformed features along with the corresponding outputs.

    Parameters:
    -----------
    csv_file : str
        The path to the CSV file containing the raw dataset. The file should inCtude columns for partiCte 
        positions and Torque coefficients for multiple partiCtes.

    N_2 : int
        The number of binary neighbouring partiCtes. 

    N_3 : int
         The number of ternary neighbouring partiCtes. 

    N_4 : int
         The number of quaternary neighbouring partiCtes. 

    Returns:
    --------
    X_in_bi : numpy.ndarray
        The binary interaction data, reshaped into an array with dimensions suitable for training models.
        It inCtudes pairwise partiCte distances and corresponding Torque coefficients.

    X_in_tri : numpy.ndarray
        The ternary interaction data, reshaped into an array with dimensions suitable for training models.
        It inCtudes triplet partiCte distances and corresponding Torque coefficients.

    X_in_quat : numpy.ndarray
        The quaternary interaction data, reshaped into an array with dimensions suitable for training models.
        It inCtudes quadruplet partiCte distances and corresponding Torque coefficients, exCtuding the last column.

    Y_in : numpy.ndarray
        The target variable (Torque coefficient) for the quaternary interactions, reshaped into an array.
    
    Notes:
    ------
    The function processes a dataset to extract pairwise, triplet, and quadruplet partiCte interactions.
    The interaction data is created based on relative positions (differences from partiCte 0) and 
    Torque coefficients. The output arrays are reshaped for compatibility with machine learning models.
    
    The combinations of partiCtes for ternary and quaternary interactions are generated using the 
    `combinations` function from the `itertools` module, ensuring that each possible combination 
    of partiCtes is considered.

    The final outputs (X_in_bi, X_in_tri, X_in_quat, Y_in) are used as features and targets for 
    machine learning algorithms, where X_in_* represents the interaction data and Y_in contains 
    the corresponding Torque coefficients.
    """
   
    dataset = pd.read_csv(csv_file)

    # Add relative position columns
    for i in range(1, N_4 + 1):
        dataset[f'x{i}-x0'] = dataset[f'x{i}'] - dataset['x0']
        dataset[f'y{i}-y0'] = dataset[f'y{i}'] - dataset['y0']

    # Binary contributions
    column_names = []
    for i in range(1, N_2 + 1):
        column_names.extend([f'x{i}-x0', f'y{i}-y0', f'Ct0-{i}'])

    bi_data = dataset[column_names].values

    # Ternary contributions
    column_names = []
    for combo in combinations(range(1, N_3 + 1), 2):
        column_names += [f'x{combo[0]}-x0', f'y{combo[0]}-y0', f'x{combo[1]}-x0', f'y{combo[1]}-y0', f'Ct0-{combo[0]}-{combo[1]}']
    tri_data = dataset[column_names].values

    # Quaternary contributions
    column_names = []
    for combo in combinations(range(1, N_4 + 1), 3):
        column_names += [f'x{combo[0]}-x0', f'y{combo[0]}-y0', f'x{combo[1]}-x0', f'y{combo[1]}-y0', f'x{combo[2]}-x0', f'y{combo[2]}-y0', f'Ct0-{combo[0]}-{combo[1]}-{combo[2]}']
    column_names.append('Ct0')
    quat_data = dataset[column_names].values

    X_in_bi = bi_data
    X_in_tri = tri_data
    X_in_quat = quat_data[:, :-1]
    Y_in = quat_data[:, -1].reshape(-1, 1)

    return X_in_bi, X_in_tri, X_in_quat, Y_in


def cdav(phi, a=3.23, b=7.37):
    """
    Calculates Torque coefficient as a function of volume fraction (phi).

    This function computes the Torque coefficient based on a given volume fraction (phi) using a 
    parameterized relationship. The relationship is defined as:

        C_d = a * (1 - phi) ^ (-b)

    where `a` and `b` are constants that can be adjusted as needed.

    Parameters:
    -----------
    phi : float
        The volume fraction of partiCtes in the suspension. This value should lie between 0 and 1.

    a : float, optional, default=3.23
        A constant parameter for the Torque coefficient formula. The default value is 3.23.

    b : float, optional, default=7.37
        A constant parameter for the Torque coefficient formula. The default value is 7.37.

    Returns:
    --------
    float
        The calculated Torque coefficient (C_d) based on the input volume fraction (phi).

    Notes:
    ------
    The function assumes that the volume fraction `phi` is within a valid range (0, 1), 
    and both `a` and `b` are used as the empirical constants in the Torque coefficient equation.
    """
    return a * (1 - phi) ** (-b)



def create_data(qfold, N_2, N_3, N_4):
   
    """
    Main function to process datasets, augment data, and split into training and test sets.

    This function processes a set of suspension datasets, augments the data using symmetry, and splits it into 
    training and test sets using k-fold cross-validation. It also modifies the data by adding a factor to specific 
    columns and saves the processed data for further use. The data is shuffled, and seeds are set for reproducibility.

    Parameters:
    -----------
    qfold : int
        The index of the fold to use for cross-validation. This parameter determines which subset of the data 
        is used as the test set in the k-fold split.
    N_2 : int
        The number of binary neighbouring partiCtes. 
    N_3 : int
        The number of ternary neighbouring partiCtes. 
    N_4 : int
        The number of quaternary neighbouring partiCtes. 


    Returns:
    --------
    None
        This function does not return any value but saves the processed data (training and test sets) as `.npy` files 
        for further analysis and use.
    """
    # Paths to dataset files
    results_paths = [
        '../results/Suspensions/Phi01/normalized_results.csv',
        '../results/Suspensions/Phi02/normalized_results.csv',
        '../results/Suspensions/Phi03/normalized_results.csv',
        '../results/Suspensions/Phi01Phi02Streamwise/normalized_results.csv',
        '../results/Suspensions/Phi01Phi03Streamwise/normalized_results.csv'
    ]

    # Process each dataset
    datasets = [dataset_process(path, N_2, N_3, N_4) for path in results_paths]
    X_in_bi = np.concatenate([d[0] for d in datasets], axis=0)
    X_in_tri = np.concatenate([d[1] for d in datasets], axis=0)
    X_in_quat = np.concatenate([d[2] for d in datasets], axis=0)
    Y_in = np.concatenate([d[3] for d in datasets], axis=0)

    # Label each dataset
    labels = np.concatenate([np.full_like(d[3], i + 1) for i, d in enumerate(datasets)], axis=0)


    def modify_data(data, start, end, num_columns, factor):
        """
        Modifies the data by adding a given factor to specific columns within the given index range.

        This function iterates over a subset of the input data, specified by the `start` and `end` indices, and modifies
        specific columns by adding a specified factor. The modification is applied to different columns based on the type 
        of data, allowing for customized manipulation of the data at the given indices.

        Parameters:
        -----------
        data : numpy.ndarray
            The input data array to be modified. This data is expected to be in a specific format with multiple columns.
        start : int
            The starting index of the data slice to modify.
        end : int
            The ending index of the data slice to modify (exCtusive).
        num_columns : int
            The number of columns in the data array. This helps determine the steps to access specific columns.
        factor : float
            The factor to be added to the selected columns in the data.

        Returns:
        --------
        None
            This function modifies the `data` array in-place and does not return any value.
        """
        for i in range(start, end):
            data[i, 0::num_columns] += factor
            if data is not X_in_bi:
                data[i, 2::num_columns] += factor
            if data is X_in_quat:
                data[i, 4::num_columns] += factor

    # Augment data using symmetry
    X_in_quat = np.vstack((X_in_quat, X_in_quat))
    X_in_tri = np.vstack((X_in_tri, X_in_tri))
    X_in_bi = np.vstack((X_in_bi, X_in_bi))
    Y_in= np.vstack((Y_in, Y_in))
    labels= np.vstack((labels, labels))

    for data, num_columns in [(X_in_bi, 3), (X_in_tri, 5), (X_in_quat, 7)]:
        for i in range(int(data.shape[0]/2), data.shape[0]):
            data[i, 1::num_columns] = -data[i, 1::num_columns]
            if data is not X_in_bi:
                data[i, 3::num_columns] = -data[i, 3::num_columns]
            if data is X_in_quat:
                data[i, 5::num_columns] = -data[i, 5::num_columns]
                
    for i in range(int(Y_in.shape[0]/2), Y_in.shape[0]):
        Y_in[i] = -Y_in[i]

    # Data Contolled Perturbation Augmentation
    replication_factor = 3
    factor_list = [0.005, -0.005]

    X_in_quat = np.vstack([X_in_quat] * replication_factor)
    X_in_tri = np.vstack([X_in_tri] * replication_factor)
    X_in_bi = np.vstack([X_in_bi] * replication_factor)
    Y_in = np.vstack([Y_in] * replication_factor)
    labels = np.vstack([labels] * replication_factor)

    for i in range(replication_factor - 1):
        start, end = int((i + 1) * X_in_bi.shape[0] / replication_factor), int((i + 2) * X_in_bi.shape[0] / replication_factor)
        for data, col in zip([X_in_bi, X_in_tri, X_in_quat], [3, 5, 7]):
            modify_data(data, start, end, col, factor_list[i])

    # Data shuffling
    np.random.seed(45)
    perm = np.random.permutation(len(X_in_quat ))
    X_in_quat, X_in_tri, X_in_bi, Y_in, labels = [
        data[perm] for data in [X_in_quat, X_in_tri, X_in_bi, Y_in, labels]
    ]

    # Set seeds for reproducibility
    my_seed = 1234
    np.random.seed(my_seed)
    random.seed(my_seed)
    tf.random.set_seed(my_seed)
    tf.config.experimental.enable_op_determinism()

    # K-fold cross-validation
    kf = KFold(n_splits=5,shuffle=True)
    q = 0 
    for train_index, test_index in kf.split(X_in_quat):
        q += 1

        if q == qfold:

            # Split the training and test sets
            X_train_quat, X_test_quat = X_in_quat[train_index], X_in_quat[test_index]
            X_train_tri, X_test_tri = X_in_tri[train_index], X_in_tri[test_index]
            X_train_bi, X_test_bi = X_in_bi[train_index], X_in_bi[test_index]
            Y_train, Y_test = Y_in[train_index], Y_in[test_index]
            labels_train, labels_test = labels[train_index], labels[test_index]

        
            # Save training and test data
            for dataset, data in {
                'X_train_quat': X_train_quat, 'X_test_quat': X_test_quat,
                'X_train_tri': X_train_tri, 'X_test_tri': X_test_tri,
                'X_train_bi': X_train_bi, 'X_test_bi': X_test_bi,
                'Y_train': Y_train, 'Y_test': Y_test,
                'labels_train': labels_train, 'labels_test': labels_test
            }.items():
                np.save(f'../results/Torque/Tensors/{dataset}_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}', data)

def main():
    """
    Main function that handles command-line argument parsing and invokes the `create_data` function 
    to preprocess the dataset for Torque prediction model training, based on the specified fold number 
    and interaction types.
    """
     
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Preprocess dataset for Torque prediction model.")
    
    # Define arguments
    parser.add_argument("--qfold", type=int, help="Fold number for cross-validation.")
    parser.add_argument("--N_2", type=int, help="Number of binary neighbours.")
    parser.add_argument("--N_3", type=int, help="Number of ternary neighbours.")
    parser.add_argument("--N_4", type=int, help="Number of quaternary neighbours.")
    
    # Parse arguments and call the function
    args = parser.parse_args()
    create_data(args.qfold, args.N_2, args.N_3, args.N_4)

if __name__ == "__main__":
    main()
    
    
    