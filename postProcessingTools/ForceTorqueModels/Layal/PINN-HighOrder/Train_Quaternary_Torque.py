# A Physics-Informed Neural Network for Modeling Higher Order Hydrodynamic Interactions in Heterogeneous Suspensions
# University of British Columbia, August 2024.
# Layal Jbara
#
# This file, Train_Quaternary_Torque.py, contains functions for training a Quaternary regression Torque model.
# The script handles loading the  preprocessed training and testing datasets.It then builds, compiles, 
# and trains a neural network model.

# Key functionalities include:
# - Model definition and layer customization
# - Implementation of callbacks for saving the best model, early stopping, learning rate decay, and tensorboard logging
# - Model evaluation on a separate test set after training
# - Handling multiple training folds for cross-validation
# - Outputting training logs and results to track the model's performance over time.

#####################################################################################################
"""
Train_Quaternary_Torque.py contains functions for training a Quaternary regression Torque model.
"""
#####################################################################################################

import os
import random
import time
import warnings
import numpy as np
import tensorflow as tf
from itertools import combinations
from math import comb
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from keras import metrics
from tensorflow.keras import initializers, backend as K
from tensorflow.keras.callbacks import Callback, EarlyStopping, ReduceLROnPlateau, TensorBoard, ModelCheckpoint
from tensorflow.keras.layers import Concatenate, Dense, Dropout, Input, Lambda, BatchNormalization, Add
from tensorflow.keras.models import Model
from tensorflow.keras.utils import set_random_seed
from keras.models import load_model
import argparse

# Set logging level and filter warnings
tf.get_logger().setLevel('ERROR')
warnings.filterwarnings('ignore')

# Print the number of available GPUs
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

try:
    # Disable all GPUs
    tf.config.set_visible_devices([], 'GPU')
except Exception as e:
    print(f"Error occurred: {e}")

#####################################################################################################

def data_train(qfold, N_2, N_3, N_4, B_S, Cds=2.74931):
    # Number of neighbors for each interaction level
    num_nghb_2 = comb(N_2, 1)  # Binary neighbors
    num_nghb_3 = comb(N_3, 2)  # Ternary neighbors
    num_nghb_4 = comb(N_4, 3)  # Quaternary neighbors

    # Hyper-parameters
    N_neurons_quat = 40
    N_neurons_tri = 30
    N_neurons_bi = 20
    epochs = 1000
    lr = 5e-4
    batch_size = B_S
    verbose = 2
    
    # Callbacks: logging and Checkpoint
    log_dir = f"../results/logs/Torque/Quaternary_q_{qfold}_N2_{N_2}_N3_{N_3}_N4_{N_4}_BS_{B_S}"
    log_path = f"../results/logs/Torque/Quaternary_log_q_{qfold}_N2_{N_2}_N3_{N_3}_N4_{N_4}_BS_{B_S}.txt"
    cp_filepath = f'../results/Checkpoints/Torque/CP-Quaternary_model_fold_{qfold}_N2_{N_2}_N3_{N_3}_N4_{N_4}_BS_{B_S}.h5'
    tensorboard_callback = TensorBoard(log_dir=log_dir, histogram_freq=1)
    Checkpoint = ModelCheckpoint(
        cp_filepath, monitor='val_loss', verbose=1,
        save_best_only=True, mode='min', period=100
    )

    # Remove existing log file if it exists
    if os.path.exists(log_path):
        os.remove(log_path)

    # Custom callback to log training results every 10 epochs
    class PrintTrainingOnTextEvery10EpochsCallback(Callback):
        def __init__(self, logpath):
            self.logpath = logpath

        def on_epoch_end(self, epoch, logs=None):
            """Logs training metrics to a text file every 10 epochs."""
            with open(self.logpath, "a") as writefile:
                if epoch == 0:
                    writefile.write("\n")
                if (epoch % 10) == 0:
                    writefile.write(
                        f"Epoch: {epoch:>3} | Loss: {logs['loss']:.4e} | "
                        f"R2: {logs['R2']:.4e} | Validation loss: {logs['val_loss']:.4e} | "
                        f"Validation R2: {logs['val_R2']:.4e} | Learning rate: {logs['lr']:.4e}\n"
                    )

    # Callback instance creation
    ten_epochs = [
        PrintTrainingOnTextEvery10EpochsCallback(logpath=log_path),
    ]

    # Reproducibility setup
    my_seed = 1234
    np.random.seed(my_seed)
    random.seed(my_seed)
    tf.random.set_seed(my_seed)
    set_random_seed(my_seed)
    tf.config.experimental.enable_op_determinism()

    # K-fold cross-validation setup
    kf = KFold(n_splits=5, shuffle=True)  # 5-fold cross-validation
    q = qfold
    print(f'+++++++++++++++++++++++++ Fold {q:1d} +++++++++++++++++++++++++')


    # Load Quaternary training and testing datasets
    X_train_quat       = np.load(f'../results/Tensors/Torque/X_train_quat_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}.npy')
    X_test_quat        = np.load(f'../results/Tensors/Torque/X_test_quat_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}.npy')

    X_train_tri       = np.load(f'../results/Tensors/Torque/X_train_tri_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}.npy')
    X_test_tri        = np.load(f'../results/Tensors/Torque/X_test_tri_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}.npy')

    X_train_bi = np.load(f'../results/Tensors/Torque/X_train_bi_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}.npy')
    X_test_bi = np.load(f'../results/Tensors/Torque/X_test_bi_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}.npy')

    Y_train = np.load(f'../results/Tensors/Torque/Y_train_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}.npy')
    Y_test = np.load(f'../results/Tensors/Torque/Y_test_q_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}.npy')

    # Define initializer
    glorot = tf.keras.initializers.GlorotNormal(seed=my_seed)

    # Summation Layer
    def summation(p_tensor):
        return tf.reduce_sum(p_tensor, axis=1, keepdims=True)
    #  ===================== Hidden layers - Quaternary================================
    quat_layer_1 = Dense(
        N_neurons_quat,
        activation='tanh',
        kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01),
        name='quat_layer_1'
    )

    quat_layer_2 = Dense(
        N_neurons_quat,
        activation='tanh',
        kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01),
        name='quat_layer_2'
    )

    quat_layer_3 = Dense(
        N_neurons_quat,
        activation='tanh',
        kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01),
        name='quat_layer_3'
    )
    
    # ===================== Processing Quaternary Neighbors with unified functionality =====================
    X_train_list_quat=[]
    input_list_quat = []
    out_list_quat = []
    for p in range(num_nghb_4):
        X_train_list_quat.append(X_train_quat[:, 7 * p:7 * p + 7])
        input_list_quat.append(Input(shape=7, name=f'quat_nghb_{p}'))
        out_list_quat.append(quat_layer_3(quat_layer_2(quat_layer_1(input_list_quat[p]))))

        
    #  ===================== Hidden layers - Ternary================================
    tri_layer_1 = Dense(
        N_neurons_tri,
        activation='tanh',
        kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01),
        name='tri_layer_1'
    )

    tri_layer_2 = Dense(
        N_neurons_tri,
        activation='tanh',
        kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01),
        name='tri_layer_2'
    )

    tri_layer_3 = Dense(
        N_neurons_tri,
        activation='tanh',
        kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01),
        name='tri_layer_3'
    )

    # ===================== Processing Ternary Neighbors with unified functionality =====================
    X_train_list_tri=[]
    input_list_tri = []
    out_list_tri = []
    for p in range(num_nghb_3):
        X_train_list_tri.append(X_train_tri[:, 5 * p:5 * p + 5])
        input_list_tri.append(Input(shape=5, name=f'tri_nghb_{p}'))
        out_list_tri.append(tri_layer_3(tri_layer_2(tri_layer_1(input_list_tri[p]))))

        
    #  ==================== Hidden layers - Binary =====================
    bi_layer_1 = Dense(
        N_neurons_bi, activation='tanh', kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01), name='bi_layer_1'
    )
    bi_layer_2 = Dense(
        N_neurons_bi, activation='tanh', kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01), name='bi_layer_2'
    )
    bi_layer_3 = Dense(
        N_neurons_bi, activation='tanh', kernel_initializer=glorot,
        bias_initializer=initializers.Constant(0.01), name='bi_layer_3'
    )

    # ===================== Processing Binary Neighbors with unified functionality =====================
    X_train_list_bi=[]
    input_list_bi = []
    out_list_bi = []
    for p in range(num_nghb_2):
        X_train_list_bi.append(X_train_bi[:, 3 * p : 3 * p + 3])
        input_list_bi.append(Input(shape=3, name=f'bi_nghb_{p}'))
        out_list_bi.append(bi_layer_3(bi_layer_2(bi_layer_1(input_list_bi[p]))))

    # Concatenate outputs from all neighbors
    concat = Concatenate()([Lambda(summation, name=f'quat_summation{j}')(out_list_quat[j]) for j in range(num_nghb_4)]+
                             [Lambda(summation, name=f'tri_summation{j}')(out_list_tri[j]) for j in range(num_nghb_3)]+
                             [Lambda(summation, name=f'bi_summation{j}')(out_list_bi[j]) for j in range(num_nghb_2)])


    # Output Layer
    nonneg = tf.keras.constraints.NonNeg()
    output_layer = Dense(1, activation='linear', kernel_constraint=nonneg, name='out_x', use_bias=False)(concat)

    # Evaluation Metric: R^2 Score
    def R2(y_true, y_pred):
        SS_res = K.sum(K.square(y_true - y_pred))
        SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
        return 1 - SS_res / (SS_tot + K.epsilon())

    # Neural Network Definition
    model = Model(inputs=[input_list_quat]+[input_list_tri]+[input_list_bi], outputs=output_layer)

    # Optimizer
    myoptimizer = tf.keras.optimizers.legacy.Adam(learning_rate=lr)

    # Early Stopping and Learning Rate Decay
    EarlyStop = EarlyStopping(monitor='val_loss', patience=50, min_delta=1e-5, start_from_epoch=200, restore_best_weights=True)
    LRDecay = ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=10, verbose=verbose)

    # Model Compilation
    model.compile(optimizer=myoptimizer, loss='mse', metrics=[R2])

    # Model Fitting
    if os.path.exists(cp_filepath):
        # Load pre-trained model if checkpoint exists
        new_model = load_model(cp_filepath, custom_objects={'R2': R2})
        history = new_model.fit(
            [X_train_list_quat]+[X_train_list_tri]+[X_train_list_bi], Y_train,
            validation_split=0.2, shuffle=True, epochs=epochs,
            verbose=verbose, batch_size=batch_size,
            callbacks=[EarlyStop, LRDecay, ten_epochs, Checkpoint, tensorboard_callback],
        )
    else:
        # Train the model from scratch
        history = model.fit(
            [X_train_list_quat]+[X_train_list_tri]+[X_train_list_bi], Y_train,
            validation_split=0.2, shuffle=True, epochs=epochs,
            verbose=verbose, batch_size=batch_size,
            callbacks=[EarlyStop, LRDecay, ten_epochs, Checkpoint, tensorboard_callback],
        )
        
    # Saving Model
    model.save(f"../results/Models/Torque/Quaternary_model_fold_{q}_N2_{N_2}_N3_{N_3}_N4_{N_4}_BS_{B_S}.h5")

    # ========================== Training Set Results ==========================
    start_time = time.perf_counter()
    Y_pred = model.predict([X_train_list_quat]+[X_train_list_tri]+[X_train_list_bi])
    end_time = time.perf_counter()
    train_R2 = r2_score(Y_train, Y_pred)
    train_time=end_time - start_time
    print(f'\033[1m\x1b[31mTraining R^2\t= {train_R2:.3f}\x1b[0m\033[0m')
    print(f'CPU time: {train_time:.3f} seconds for {Y_pred.shape[0]:8d} datapoints')

    # =========================== Test Set Results ============================
    # Prepare the test set
    X_test_list_bi=[]
    X_test_list_tri=[]
    X_test_list_quat=[]
    X_test_list_bi = [X_test_bi[:, 3 * p:3 * p + 3] for p in range(num_nghb_2)]
    X_test_list_tri= [X_test_tri[:, 5 * p:5 * p + 5] for p in range(num_nghb_3)]
    X_test_list_quat= [X_test_quat[:, 7 * p:7 * p + 7] for p in range(num_nghb_4)]
    start_time = time.perf_counter()
    Y_pred_test = model.predict([X_test_list_quat]+[X_test_list_tri]+[X_test_list_bi])
    end_time = time.perf_counter()
    test_R2 = r2_score(Y_test, Y_pred_test)
    test_time=end_time - start_time
    print(f'\033[1m\x1b[31mTest R^2\t= {test_R2:.3f}\x1b[0m\033[0m')
    print(f'CPU time: {test_time:.3f} seconds for {Y_pred_test.shape[0]:8d} datapoints')

    # ============================= Summary ===================================
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('========================================\nOverall Performance:')
    print(f'<Train R^2> = {train_R2:.3f}')
    print(f'<Test R^2> = {test_R2:.3f}')
    print(f'<Train Time> = {train_time:.3f} seconds')
    print(f'<Test Time> = {test_time:.3f} seconds')

    # Model details
    print(f'\nTotal number of model parameters: {model.count_params()}\n')
    print(model.summary())
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qfold', type=int, required=True, help='Fold number for cross-validation')
    parser.add_argument('--N_2', type=int, required=True, help='Number of Binary neighbours')
    parser.add_argument('--N_3', type=int, required=True, help='Number of Ternary neighbours')
    parser.add_argument('--N_4', type=int, required=True, help='Number of Quaternary neighbours')
    parser.add_argument('--B_S', type=int, required=True, help='Batch size for training')
    
    # Parse arguments and call the traning function
    args = parser.parse_args()
    data_train(args.qfold, args.N_2, args.N_3, args.N_4, args.B_S)
    print(f"Training completed. Model saved!")
    
if __name__ == "__main__":
    main()