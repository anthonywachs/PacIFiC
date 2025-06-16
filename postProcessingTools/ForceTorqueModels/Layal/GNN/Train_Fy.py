"""
Graph Neural Networks for Modeling Hydrodynamic Closure Laws in Non-Spherical Particle-Laden Flows
University of British Columbia, March 2025
Author: Layal Jbara

This script, Train_Fy.py, includes functions that:
- Convert data into a graph structure.
- Design and train Graph Neural Networks (GNNs) for a regression-based task to predict the hydrodynamic lift forces.
"""

#####################################################################################################
# Import Libraries
#####################################################################################################

import os
import json
import random
import argparse
import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial.distance import cdist
from sklearn.metrics import r2_score

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch_geometric
import torch_scatter
from torch.optim.lr_scheduler import StepLR
from torch_scatter import scatter_add
from torch_geometric.nn import MessagePassing, global_mean_pool
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.utils import softmax
import training_functions as trf
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

#####################################################################################################
# Physics-Informed Graph Neural Network (GNN) Layer
#####################################################################################################

class PhysicsInformedGNN(MessagePassing):
    """
    A physics-informed Graph Neural Network (GNN) for modeling hydrodynamic closure laws.
    
    This GNN incorporates a multi-head attention mechanism.
    """
    
    def __init__(self, in_channels, edge_channels, hidden_channels, num_heads=1, activation_func=F.gelu, dropout_rate=0.4):
        super().__init__(aggr=None)  # Disable default aggregation

        self.num_heads = num_heads
        self.activation_func = activation_func

        # Dynamic attention mechanism (multi-head)
        self.attn_fc = nn.ModuleList([nn.Linear(hidden_channels, 1) for _ in range(num_heads)])

        # Message passing layers
        self.msg1 = nn.Linear(2 * in_channels + edge_channels, hidden_channels)
        self.msg2 = nn.Linear(hidden_channels, hidden_channels)
        self.msg3 = nn.Linear(hidden_channels, hidden_channels)
        self.msg4 = nn.Linear(hidden_channels, hidden_channels)

        # Dropout layers to prevent overfitting
        self.dropout1 = nn.Dropout(p=dropout_rate)
        self.dropout2 = nn.Dropout(p=dropout_rate)
        self.dropout3 = nn.Dropout(p=dropout_rate)
        self.dropout4 = nn.Dropout(p=dropout_rate)

        # Initialize weights
        self.reset_parameters()

    def reset_parameters(self):
        """Initialize model weights using Xavier uniform distribution."""
        nn.init.xavier_uniform_(self.msg1.weight)
        nn.init.xavier_uniform_(self.msg2.weight)
        nn.init.xavier_uniform_(self.msg3.weight)
        nn.init.xavier_uniform_(self.msg4.weight)
        for attn_layer in self.attn_fc:
            nn.init.xavier_uniform_(attn_layer.weight)

    def forward(self, x, edge_index, edge_attr):
        """Forward pass propagating node features through the graph."""
        return self.propagate(edge_index, x=x, edge_attr=edge_attr)

    def message(self, x_i, x_j, edge_attr):
        """Compute message passing updates for node pairs."""
        # Concatenate node and edge features
        msg_input = torch.cat([x_i, x_j, edge_attr], dim=1)

        # Apply multiple message layers with activation
        msg_out = self.activation_func(self.msg1(msg_input))
        msg_out = self.activation_func(self.msg2(msg_out))
        msg_out = self.activation_func(self.msg3(msg_out))

        return msg_out
    
    def aggregate(self, inputs, index, dim_size):
        """Aggregate incoming messages using multi-head attention."""
        # Compute attention scores for each head
        attn_scores = [softmax(F.leaky_relu(attn_layer(inputs), negative_slope=0.1), index=index, num_nodes=dim_size) 
                       for attn_layer in self.attn_fc]
        
        # Average attention weights across heads
        attn_weights = torch.mean(torch.stack(attn_scores), dim=0)

        # Perform weighted aggregation
        aggr_out = scatter_add(attn_weights * inputs, index, dim=0, dim_size=dim_size)

        return aggr_out

    def update(self, aggr_out):
        """Apply activation function after aggregation."""
        return self.activation_func(aggr_out)

#####################################################################################################
# Sequential Graph Neural Network (GNN) Class
#####################################################################################################

class SequentialPhysicsInformedGNN(nn.Module):
    """
    A sequential Graph Neural Network (GNN) model with physics-informed layers.
    This model applies multiple GNN layers in sequence, followed by fully connected layers.
    """
    def __init__(self, in_channels, edge_channels, hidden_channels, out_channels, num_heads=1, activation_func=F.gelu, dropout_rate=0.5):
        super().__init__()

        # Define multiple GNN layers applied sequentially
        self.gnn1 = PhysicsInformedGNN(in_channels, edge_channels, hidden_channels, num_heads, activation_func)
        self.gnn2 = PhysicsInformedGNN(hidden_channels, edge_channels, hidden_channels, num_heads, activation_func)
        self.gnn3 = PhysicsInformedGNN(hidden_channels, edge_channels, hidden_channels, num_heads, activation_func)
        self.gnn4 = PhysicsInformedGNN(hidden_channels, edge_channels, hidden_channels, num_heads, activation_func)

        # Fully connected (FC) layers for final prediction
        self.fc1 = nn.Linear(hidden_channels, hidden_channels // 2)
        self.fc2 = nn.Linear(hidden_channels // 2, hidden_channels // 4)
        self.fc3 = nn.Linear(hidden_channels // 4, out_channels)

        # Activation function and dropout layers
        self.activation_func = activation_func
        self.dropout = nn.Dropout(p=dropout_rate)
        self.dropout2 = nn.Dropout(p=dropout_rate)

        # Initialize weights of the fully connected layers
        self.reset_parameters()

    def reset_parameters(self):
        """Initialize fully connected layer weights using Xavier uniform initialization."""
        nn.init.xavier_uniform_(self.fc1.weight)
        nn.init.xavier_uniform_(self.fc2.weight)

    def forward(self, x, edge_index, edge_attr):
        """
        Forward pass of the sequential physics-informed GNN model.
        
        Args:
            x (Tensor): Node feature matrix.
            edge_index (Tensor): Graph connectivity (edge list).
            edge_attr (Tensor): Edge feature matrix.
        
        Returns:
            Tensor: Output predictions.
        """
        # Apply GNN layers sequentially
        x = self.gnn1(x, edge_index, edge_attr)
        x = self.gnn2(x, edge_index, edge_attr)
        x = self.gnn3(x, edge_index, edge_attr)
        x = self.gnn4(x, edge_index, edge_attr)

        # Fully connected layers with activation and dropout
        x = self.fc1(x)
        x = self.activation_func(x)
        x = self.dropout(x)

        x = self.fc2(x)
        x = self.activation_func(x)
        x = self.dropout2(x)

        x = self.fc3(x)  # Output layer
        return x



#####################################################################################################
# GNN Training Function 
#####################################################################################################

def train_model_four_sequential(
    directory_name,
    hidden_channels,
    dataloader,
    N_neigh,
    num_epochs=100,
    learning_rate=0.01,
    patience=50,
    verbose=True
):
    """
    Trains a sequential physics-informed GNN model with early stopping.

    Args:
        directory_name (str): Directory to save results.
        hidden_channels (int): Number of hidden layer channels.
        dataloader: DataLoader for training and validation.
        N_neigh (int): Number of neighbors.
        num_epochs (int, optional): Maximum number of training epochs. Defaults to 100.
        learning_rate (float, optional): Learning rate for the optimizer. Defaults to 0.01.
        patience (int, optional): Number of consecutive increases in validation loss allowed before stopping. Defaults to 5.
        verbose (bool, optional): If True, prints progress updates. Defaults to True.

    Returns:
        dict: Best model state dictionary.
    """
    # Set device to GPU if available, otherwise use CPU
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Initialize results storage
    results = []
    json_filename = f"{directory_name}Fy_det_results_{hidden_channels}channels_neighbours_{N_neigh}_optm.json"
    
    # Load existing results if the JSON file exists
    if os.path.exists(json_filename):
        with open(json_filename, "r") as f:
            results = json.load(f)

    # Initialize model and move to appropriate device
    model = SequentialPhysicsInformedGNN(
        in_channels=2,
        edge_channels=3,
        hidden_channels=hidden_channels,
        out_channels=1,
        num_heads=5,
        dropout_rate=0.2
    ).to(device)
    model = model.double()

    # Define optimizer with weight decay
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=1e-5)

    # Initialize metrics storage
    loss_history, valid_loss_history = [], []
    r2_train_history, r2_val_history = [], []

    # Early stopping variables
    best_valid_loss = 0
    best_model_state = None
    consec_increases = 0

    for epoch in range(num_epochs):
        total_loss, valid_total_loss = 0, 0
        
        all_targets, all_predictions = [], []
        train_targets, train_predictions = [], []
        valid_targets, valid_predictions = [], []
        
        for data in dataloader:
            model.train()
            
            # Move data to device
            data.x = data.x.to(torch.float64).to(device)
            data.edge_index = data.edge_index.to(device)
            data.edge_attr = data.edge_attr.to(torch.float64).to(device)
            data.y = data.y.to(torch.float64).to(device)
            
            # Forward pass and loss computation
            optimizer.zero_grad()
            out = model(data.x, data.edge_index, data.edge_attr)
            train_loss = F.mse_loss(out[data.train_mask], data.y[data.train_mask])
            train_loss.backward()
            optimizer.step()
            total_loss += train_loss.item()
            
            # Validation phase
            model.eval()
            with torch.no_grad():
                valid_out = model(data.x, data.edge_index, data.edge_attr)
                valid_loss = F.mse_loss(valid_out[data.val_mask], data.y[data.val_mask])
                valid_total_loss += valid_loss.item()
                
                # Store targets and predictions
                all_targets.append(data.y.cpu().numpy())
                train_targets.append(data.y[data.train_mask].cpu().numpy())
                valid_targets.append(data.y[data.val_mask].cpu().numpy())
                
                all_predictions.append(out.cpu().numpy())
                train_predictions.append(out[data.train_mask].cpu().numpy())
                valid_predictions.append(out[data.val_mask].cpu().numpy())

        # Compute epoch metrics
        avg_train_loss = total_loss / len(dataloader)
        avg_valid_loss = valid_total_loss / len(dataloader)
        loss_history.append(avg_train_loss)
        valid_loss_history.append(avg_valid_loss)

        # Concatenate stored data
        all_targets, train_targets, valid_targets = map(np.concatenate, [all_targets, train_targets, valid_targets])
        all_predictions, train_predictions, valid_predictions = map(np.concatenate, [all_predictions, train_predictions, valid_predictions])
        
        # Compute R² scores
        r2_tot = r2_score(all_targets, all_predictions)
        r2_train = r2_score(train_targets, train_predictions)
        r2_valid = r2_score(valid_targets, valid_predictions)
        
        if verbose and epoch % 50 == 0:
            print(f"Epoch {epoch + 1}/{num_epochs}, Train Loss: {avg_train_loss:.4f}, Valid Loss: {avg_valid_loss:.4f}, "
                  f"R² Train: {r2_train:.4f}, R² Valid: {r2_valid:.4f}, lr: {optimizer.param_groups[0]['lr']:.4f}")
        
        # Early stopping logic
        if r2_valid > best_valid_loss:
            best_valid_loss = r2_valid
            best_model_state = model.state_dict()
            consec_increases = 0
            best_r2_train, best_r2_valid, best_r2_total = r2_train, r2_valid, r2_tot
            best_model, best_epoch = model, epoch
        else:
            consec_increases += 1
        
        if consec_increases >= patience and epoch > 2000:
            if verbose:
                print(f"Early stopping triggered at epoch {epoch + 1}")
            break
    
    # Store configuration results
    config_result = {
        "hidden_channels": hidden_channels,
        "Neighbours": N_neigh,
        "num_epochs": epoch,
        "final_r2_train": r2_train,
        "final_r2_valid": r2_valid,
        "final_r2_total": r2_tot,
        "best_epoch": best_epoch,
        "best_r2_train": best_r2_train,
        "best_r2_valid": best_r2_valid,
        "best_r2_total": best_r2_total
    }
    results.append(config_result)
    
    # Save results to JSON
    try:
        with open(json_filename, "w") as f:
            json.dump(results, f)
            print("Results saved to JSON file.")
    except Exception as e:
        print("Failed to write JSON file:", e)
    
    # Save the best model state to a file
    model_save_path = f"{directory_name}best_model_{hidden_channels}channels_neighbours_{N_neigh}_optm.pt"
    torch.save(best_model_state, model_save_path)
    print(f"Best model saved to {model_save_path}")

    return best_model_state



def parse_args():
    """
    Parse command-line arguments for configuring the neural network.

    Returns:
    argparse.Namespace: Parsed arguments containing model parameters.
    """
    parser = argparse.ArgumentParser(description="Model parameters for neural network.")

    # Model hyperparameters
    parser.add_argument('--N_neigh', type=int, default=10, 
                        help="Number of neighbors. Default is 10.")
    parser.add_argument('--directory_name', type=str, required=True, 
                        help="Directory name for saving model outputs.")
    parser.add_argument('--hidden_channels', type=int, default=64, 
                        help="Number of hidden channels in the network. Default is 64.")

    return parser.parse_args()


if __name__ == "__main__":

    # Parse command-line arguments
    args = parse_args()

    # Display parsed arguments
    print("Number of Neighbors:", args.N_neigh)
    print("Directory Name:", args.directory_name)
    print("Hidden Channels:", args.hidden_channels)

    # Set the random seed for reproducibility
    trf.set_seed(34695)

    # Set device for computation (GPU if available, otherwise CPU)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Load CSV files for different dataset variations
    df = pd.read_csv(args.directory_name + 'boundary_neighbours_new_scale.csv')
    df2 = pd.read_csv(args.directory_name + 'angle.csv')

    dfy = pd.read_csv(args.directory_name + 'boundary_neighbours_new_scale_flip_y.csv')
    df2y = pd.read_csv(args.directory_name + 'angle_flip_y.csv')

    dfz = pd.read_csv(args.directory_name + 'boundary_neighbours_new_scale_flip_z.csv')
    df2z = pd.read_csv(args.directory_name + 'angle_flip_z.csv')

    dfzy = pd.read_csv(args.directory_name + 'boundary_neighbours_new_scale_flip_zy.csv')
    df2zy = pd.read_csv(args.directory_name + 'angle_flip_yz.csv')



        # Function to prepare dataset for a given flipped configuration
    def prepare_data(df, df2, force='Cl0', flip_y=False):
        """
        Prepares the dataset by applying transformations, optionally flipping the force component,
        and moving it to the specified device.

        Parameters:
        ----------
        df : pd.DataFrame
            Data containing boundary neighbors.
        df2 : pd.DataFrame
            Data containing angles.
        force : str, optional
            Name of the force component to be used (default is 'Cl0').
        flip_y : bool, optional
            If True, the selected force component will be augumented using symmetry
            before creating the dataset. Default is False.

        Returns:
        -------
        torch_geometric.data.Data
            Processed dataset ready for training.
        """
        # Make a copy to avoid modifying the original dataframe
        new_df = df.copy()

        # Optionally flip the specified force component
        if flip_y:
            new_df[force] = -new_df[force]

        # Create and preprocess dataset
        data = trf.create_data(new_df, df2, N=args.N_neigh, RE=args.RE, force=force).to(device)
        data = trf.create_masks(data, train_ratio=0.80, seed=7264).to(device)

        return data

    # Prepare datasets for original and flipped configurations
    data_original = prepare_data(df, df2)
    data_y_flipped = prepare_data(dfy, df2y, flip_y=True)
    data_z_flipped = prepare_data(dfz, df2z)
    data_zy_flipped = prepare_data(dfzy, df2zy, flip_y=True)

    # Combine all datasets into a list for training
    all_data = [data_original, data_y_flipped, data_z_flipped, data_zy_flipped]

    # Define DataLoader with the custom collate function
    dataloader = DataLoader(all_data, batch_size=1, shuffle=True, collate_fn=trf.collate_fn)

    # Train the model using the specified arguments
    model_best_state = train_model_four_sequential(
        directory_name=args.directory_name,
        hidden_channels=args.hidden_channels,
        dataloader=dataloader,
        N_neigh=args.N_neigh,
        num_epochs=2000,
        learning_rate=0.001,
        patience=50,
        verbose=True
    )
