# A Physics-Informed Neural Network for Modeling Higher Order Hydrodynamic Interactions in Heterogeneous Suspensions
# University of British Columbia, August 2024.
# Layal Jbara
#
# This file, Inference.py, contains functions that leverage pre-trained decoupled models 
# to compute interaction terms specific to the local microstructure of each particle in the suspension.
# The primary goal of these computations is to evaluate binary, ternary, and quaternary hydrodynamic interactions
# within a suspension.

#####################################################################################################

"""
Inference.py contains functions that utilize pre-trained decoupled models to compute the interaction terms 
specific to the local microstructure of each particle in a suspension. 
"""

#####################################################################################################


import numpy as np
import pandas as pd
import timeit
import warnings
from math import comb
import multiprocessing as mp
from functools import partial
import tensorflow as tf
from keras.models import load_model
from tensorflow.keras import backend as K
from itertools import combinations
import ray
import glob
import os
import sys

# Add custom path for module imports
sys.path.append('/zfs/users/lmj44/lmj44/anaconda3/lib/python3.8/site-packages')

# Disable GPU usage (if necessary)
try:
    # Disable all GPUs
    tf.config.set_visible_devices([], 'GPU')
    visible_devices = tf.config.get_visible_devices()
    for device in visible_devices:
        assert device.device_type != 'GPU'
except:
    # Invalid device or cannot modify virtual devices once initialized
    pass


def load_models(hydro, fold_number):
    """
    This function loads pre-trained models for different interaction types (Binary, Ternary, Quaternary)
    for a given fold number in a specific hydrodynamic configuration.

    Parameters:
    - hydro: str. The hydrodynamic configuration to load models for.
    - fold_number: int. The fold number of the model to load.

    Returns:
    - models: dict. A dictionary containing the loaded models for 'Binary', 'Ternary', and 'Quaternary' interactions.
    """
    
    models = {}
    for kind in ['Binary', 'Ternary', 'Quaternary']:
        models[kind] = load_model(
            f'../../Train_Models_{hydro}/{kind}/{kind}_model_fold_{fold_number}.h5', custom_objects={'R2': R2})

    return models

def R2(y_true, y_pred):
    """
    This function computes the R² (coefficient of determination) for a set of true and predicted values.

    The R² score is a statistical measure that represents the proportion of the variance in the dependent variable 
    that is predictable from the independent variable(s). It is commonly used to evaluate the performance of regression models.

    Parameters:
    - y_true: tensor. The true values (ground truth).
    - y_pred: tensor. The predicted values from the model.

    Returns:
    - R2: tensor. The R² score, a measure of the goodness of fit, where 1 indicates perfect prediction and 0 indicates that the model is no better than using the mean of the true values.
    """
    
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    R2 = 1 - SS_res / (SS_tot + K.epsilon())
    return R2


def ray_compute_all_orders(neigh,ncores):
    """
    Initializes Ray for parallel computation and defines functions to compute particle interactions for 
    all orders of interaction (binary, ternary, quaternary). This function sets up the necessary parallel 
    processing environment and distributes the computation of pairwise (binary), triplet (ternary), and 
    quadruplet (quaternary) interactions across multiple CPU cores using Ray, improving performance for large datasets.
    
    Parameters:
    -----------
    neigh : int
        The number of neighboring particles to consider for the interaction calculations.
    ncores : int
        The number of CPU cores to be used for parallel computation with Ray.

    Returns:
    --------
    None
        This function does not return any value, as it performs computations in parallel.
    """
    
    # Get all dataset paths
    datasets = glob.iglob('../results/Suspensions/*/boundary_neighbours.csv') 
    
    for file in datasets:
        
        # Initialize Ray
        ray.shutdown()
        
        ray.init(num_cpus=ncores,log_to_driver=False)
       
        @ray.remote
        def compute_all_orders(df,interval,neigh):
            """
            This function computes interactions between particles for all orders of interaction (binary, ternary, and quaternary)
            based on the given DataFrame containing particle data.

            Parameters:
            - df (DataFrame): The input DataFrame containing particle properties such as positions and coefficients.
            - interval (list): The range of indices to process from the DataFrame.
            - neigh (int): The number of neighboring particles to consider for the computations.

            Returns:
            - df (DataFrame): The input DataFrame with additional columns for computed interaction coefficients 
                              (Drag, Lift, and Torque) for all interaction orders (binary, ternary, and quaternary).
            """

            # Select df based on interval
            df = df.iloc[interval[0]:interval[-1]]


            # Load Models
            models = load_models('Drag',1)
            models_Lift = load_models('Lift',1)
            models_Torque = load_models('Torque',1)
            
            ##_____________Binary_Computations_______________##

             # Extract Reference Particle Coordinate
            x0_values = df['x0'].values
            y0_values = df['y0'].values

            # Calculate Relative Differences 
            x_diff = df[[f'x{j}' for j in range(1, neigh + 1)]].values - x0_values[:, np.newaxis]
            y_diff = df[[f'y{j}' for j in range(1, neigh + 1)]].values - y0_values[:, np.newaxis]

            # Reshape data into input data to the model
            input_data = np.column_stack([x_diff.ravel(), y_diff.ravel()])

            # Make Binary predictions
            predictions = models['Binary'].predict(input_data, verbose=1).reshape(-1, neigh)
            predictions_Lift = models_Lift['Binary'].predict(input_data, verbose=1).reshape(-1, neigh)
            predictions_Torque = models_Torque['Binary'].predict(input_data, verbose=1).reshape(-1, neigh)

         
            # Add prediction columns to DataFrame
            for j in range(1, neigh + 1):
                df[f'Cd0-{j}'] = predictions[:, j-1]
                df[f'Cl0-{j}'] = predictions_Lift[:, j-1]
                df[f'Ct0-{j}'] = predictions_Torque[:, j-1]


            ##_____________Ternary_Computations_______________##

            # Calculate Relative Differences 
            diff_x0 = df[[f'x{combo[0]}' for combo in combinations(range(1, neigh+1), 2)]].values - x0_values[:, np.newaxis]
            diff_y0 = df[[f'y{combo[0]}' for combo in combinations(range(1, neigh+1), 2)]].values - y0_values [:, np.newaxis]
            diff_x1 = df[[f'x{combo[1]}' for combo in combinations(range(1, neigh+1), 2)]].values - x0_values[:, np.newaxis]
            diff_y1 = df[[f'y{combo[1]}' for combo in combinations(range(1, neigh+1), 2)]].values - y0_values[:, np.newaxis]

            # Extract Binary Drag Coefficients 
            cd0_combo0 = df[[f'Cd0-{combo[0]}' for combo in combinations(range(1, neigh+1), 2)]].values
            cd0_combo1 = df[[f'Cd0-{combo[1]}' for combo in combinations(range(1, neigh+1), 2)]].values
            
            # Extract Binary Lift Coefficients 
            cl0_combo0 = df[[f'Cl0-{combo[0]}' for combo in combinations(range(1, neigh+1), 2)]].values
            cl0_combo1 = df[[f'Cl0-{combo[1]}' for combo in combinations(range(1, neigh+1), 2)]].values

            # Extract Binary Torque Coefficients 
            Ct0_combo0 = df[[f'Ct0-{combo[0]}' for combo in combinations(range(1, neigh+1), 2)]].values
            Ct0_combo1 = df[[f'Ct0-{combo[1]}' for combo in combinations(range(1, neigh+1), 2)]].values


            # Reshape data into input data to the model
            input_data =[
                np.array([diff_x0.ravel(), diff_y0.ravel(), cd0_combo0.ravel()]).T,
                np.array([diff_x1.ravel(), diff_y1.ravel(), cd0_combo1.ravel()]).T
            ]

            input_data_Lift =[
                            np.array([diff_x0.ravel(), diff_y0.ravel(), cl0_combo0.ravel()]).T,
                            np.array([diff_x1.ravel(), diff_y1.ravel(), cl0_combo1.ravel()]).T
                        ]

            input_data_Torque =[
                            np.array([diff_x0.ravel(), diff_y0.ravel(), Ct0_combo0.ravel()]).T,
                            np.array([diff_x1.ravel(), diff_y1.ravel(), Ct0_combo1.ravel()]).T
                        ]

            # Make Ternary predictions
            predictions = models['Ternary'].predict(input_data, verbose=1).reshape(-1,comb(neigh, 2))
            predictions_Lift = models_Lift['Ternary'].predict(input_data_Lift, verbose=1).reshape(-1,comb(neigh, 2))
            predictions_Torque = models_Torque['Ternary'].predict(input_data_Torque, verbose=1).reshape(-1,comb(neigh, 2))

          
            # Add prediction columns to DataFrame
            for i,combo in enumerate(combinations(range(1, neigh+1), 2)):
                df[f'Cd0-{combo[0]}-{combo[1]}'] = predictions[:, i]
                df[f'Cl0-{combo[0]}-{combo[1]}'] = predictions_Lift[:, i]
                df[f'Ct0-{combo[0]}-{combo[1]}'] = predictions_Torque[:, i]
                

            ##_____________Quaternary_Computations_______________##

            # Calculate Relative Differences 
            diff_x0 = df[[f'x{combo[0]}' for combo in combinations(range(1, neigh+1), 3)]].values - x0_values[:, np.newaxis]
            diff_y0 = df[[f'y{combo[0]}' for combo in combinations(range(1, neigh+1), 3)]].values - y0_values[:, np.newaxis]
            diff_x1 = df[[f'x{combo[1]}' for combo in combinations(range(1, neigh+1), 3)]].values - x0_values[:, np.newaxis]
            diff_y1 = df[[f'y{combo[1]}' for combo in combinations(range(1, neigh+1), 3)]].values - y0_values[:, np.newaxis]
            diff_x2 = df[[f'x{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values - x0_values[:, np.newaxis]
            diff_y2 = df[[f'y{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values - y0_values[:, np.newaxis]


            # Extract Binary Drag Coefficients 
            cd0_combo0 = df[[f'Cd0-{combo[0]}' for combo in combinations(range(1, neigh+1), 3)]].values
            cd0_combo1 = df[[f'Cd0-{combo[1]}' for combo in combinations(range(1, neigh+1), 3)]].values
            cd0_combo2 = df[[f'Cd0-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values
            
            # Extract Binary Lift Coefficients 
            cl0_combo0 = df[[f'Cl0-{combo[0]}' for combo in combinations(range(1, neigh+1), 3)]].values
            cl0_combo1 = df[[f'Cl0-{combo[1]}' for combo in combinations(range(1, neigh+1), 3)]].values
            cl0_combo2 = df[[f'Cl0-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values

            # Extract Binary Torque Coefficients 
            Ct0_combo0 = df[[f'Ct0-{combo[0]}' for combo in combinations(range(1, neigh+1), 3)]].values
            Ct0_combo1 = df[[f'Ct0-{combo[1]}' for combo in combinations(range(1, neigh+1), 3)]].values
            Ct0_combo2 = df[[f'Ct0-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values

            # Extract Ternary Drag Coefficients 
            cd0_combo01 = df[[f'Cd0-{combo[0]}-{combo[1]}' for combo in combinations(range(1, neigh+1), 3)]].values
            cd0_combo02 = df[[f'Cd0-{combo[0]}-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values
            cd0_combo12 = df[[f'Cd0-{combo[1]}-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values
            
            # Extract Ternary Lift Coefficients
            cl0_combo01 = df[[f'Cl0-{combo[0]}-{combo[1]}' for combo in combinations(range(1, neigh+1), 3)]].values
            cl0_combo02 = df[[f'Cl0-{combo[0]}-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values
            cl0_combo12 = df[[f'Cl0-{combo[1]}-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values

            # Extract Ternary Torque Coefficients
            Ct0_combo01 = df[[f'Ct0-{combo[0]}-{combo[1]}' for combo in combinations(range(1, neigh+1), 3)]].values
            Ct0_combo02 = df[[f'Ct0-{combo[0]}-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values
            Ct0_combo12 = df[[f'Ct0-{combo[1]}-{combo[2]}' for combo in combinations(range(1, neigh+1), 3)]].values


            # Reshape data into input data to the model
            input_data =[   
                        np.array([diff_x0.ravel(), diff_y0.ravel(), diff_x1.ravel(), diff_y1.ravel(),cd0_combo01.ravel()]).T,
                        np.array([diff_x0.ravel(), diff_y0.ravel(), diff_x2.ravel(), diff_y2.ravel(),cd0_combo02.ravel()]).T,
                        np.array([diff_x1.ravel(), diff_y1.ravel(), diff_x2.ravel(), diff_y2.ravel(),cd0_combo12.ravel()]).T
                        ] + [
                        np.array([diff_x0.ravel(), diff_y0.ravel(),cd0_combo0.ravel()]).T,
                        np.array([diff_x1.ravel(), diff_y1.ravel(), cd0_combo1.ravel()]).T,
                        np.array([diff_x2.ravel(), diff_y2.ravel(), cd0_combo2.ravel()]).T
                        ]
            
            # Reshape data into input data to the model
            input_data_Lift =[   
                        np.array([diff_x0.ravel(), diff_y0.ravel(), diff_x1.ravel(), diff_y1.ravel(),cl0_combo01.ravel()]).T,
                        np.array([diff_x0.ravel(), diff_y0.ravel(), diff_x2.ravel(), diff_y2.ravel(),cl0_combo02.ravel()]).T,
                        np.array([diff_x1.ravel(), diff_y1.ravel(), diff_x2.ravel(), diff_y2.ravel(),cl0_combo12.ravel()]).T
                        ] + [
                        np.array([diff_x0.ravel(), diff_y0.ravel(),cl0_combo0.ravel()]).T,
                        np.array([diff_x1.ravel(), diff_y1.ravel(), cl0_combo1.ravel()]).T,
                        np.array([diff_x2.ravel(), diff_y2.ravel(), cl0_combo2.ravel()]).T
                        ]
            
            input_data_Torque =[   
                        np.array([diff_x0.ravel(), diff_y0.ravel(), diff_x1.ravel(), diff_y1.ravel(),Ct0_combo01.ravel()]).T,
                        np.array([diff_x0.ravel(), diff_y0.ravel(), diff_x2.ravel(), diff_y2.ravel(),Ct0_combo02.ravel()]).T,
                        np.array([diff_x1.ravel(), diff_y1.ravel(), diff_x2.ravel(), diff_y2.ravel(),Ct0_combo12.ravel()]).T
                        ] + [
                        np.array([diff_x0.ravel(), diff_y0.ravel(),Ct0_combo0.ravel()]).T,
                        np.array([diff_x1.ravel(), diff_y1.ravel(), Ct0_combo1.ravel()]).T,
                        np.array([diff_x2.ravel(), diff_y2.ravel(), Ct0_combo2.ravel()]).T
                        ]
            #Make Quaternary predictions
            predictions = models['Quaternary'].predict(input_data, verbose=1).reshape(-1,comb(neigh, 3))
            predictions_Lift = models_Lift['Quaternary'].predict(input_data_Lift, verbose=1).reshape(-1,comb(neigh, 3))
            predictions_Torque = models_Torque['Quaternary'].predict(input_data_Torque, verbose=1).reshape(-1,comb(neigh, 3))

            # Add prediction columns to DataFrame with a progress bar
            for i, combo in enumerate(combinations(range(1, neigh + 1), 3)):
                df[f'Cd0-{combo[0]}-{combo[1]}-{combo[2]}'] = predictions[:, i]
                df[f'Cl0-{combo[0]}-{combo[1]}-{combo[2]}'] = predictions_Lift[:, i]
                df[f'Ct0-{combo[0]}-{combo[1]}-{combo[2]}'] = predictions_Torque[:, i]

            return df

        # Load the dataset into a pandas DataFrame
        df = pd.read_csv(file)

        # Define the interval for dividing the dataset into chunks
        total_rows = len(df)
        size = total_rows // ncores
        intervals = [(i*size, (i+1)*size) for i in range(ncores)]
        if total_rows > intervals[-1][1]:
            intervals[-1] = (intervals[-1][0], total_rows)  # Add the remainder (if any)
        
        # Parallelize the computation using Ray for each chunk of data
        futures = [compute_all_orders.remote(df,interval, neigh) for interval in intervals]
        
        results = ray.get(futures)

        # Concatenate the results back into a single DataFrame
        final_df = pd.concat(results)
        final_df.to_csv(os.path.join('../results/Suspensions/',file.split(os.path.sep)[3]+'/results.csv'))
        
        # Shutdown Ray after processing
        ray.shutdown()
        
  
    
if __name__ == "__main__":  
    ray_compute_all_orders(neigh=10,ncores=36)