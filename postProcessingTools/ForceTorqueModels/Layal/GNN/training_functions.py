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


def set_seed(seed=64213):
    """
    Set the random seed for reproducibility across multiple libraries.

    Parameters:
    seed (int): The seed value for random number generators. Default is 64213.
    """
    # Set seed for Python's built-in random module
    random.seed(seed)
    
    # Set seed for NumPy
    np.random.seed(seed)
    
    # Set seed for PyTorch (CPU & GPU)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # For multi-GPU setups
    
    # Ensure deterministic behavior in PyTorch (at a potential performance cost)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)

# Define a custom collate function for batching data
def collate_fn(batch):
    """
    Custom function to stack Data objects in the batch into a single batch.
    """
    return Data(
        x=torch.cat([data.x for data in batch], dim=0),
        edge_index=torch.cat([data.edge_index for data in batch], dim=1),
        y=torch.cat([data.y for data in batch], dim=0),
        # Add other required attributes of the Data object here if needed

def get_optimizer(model, lr=1e-3, weight_decay=1e-5):
    """
    Returns an Adam optimizer with weight decay for regularization.
    
    Args:
        model (torch.nn.Module): The model to optimize.
        lr (float, optional): Learning rate. Default is 1e-3.
        weight_decay (float, optional): Weight decay for L2 regularization. Default is 1e-5.
    
    Returns:
        torch.optim.Optimizer: Adam optimizer instance.
    """
    return torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

def create_data(df, df2, N, RE, force='Cd0', device='cpu'):
    """
    Constructs a graph representation from particle data, incorporating edge features and target forces.
    
    Args:
        df (pd.DataFrame): DataFrame containing particle coordinates and their neighbors.
        df2 (pd.DataFrame): DataFrame containing rotation angles (phi, psi) for each particle.
        N (int): Number of neighbors considered per particle.
        RE (float): Reynolds number scaling factor.
        device (str, optional): Device to store tensors ('cpu' or 'cuda'). Default is 'cpu'.
    
    Returns:
        Data: PyTorch Geometric Data object containing graph information.
    """
    edge_index = []  # List to store edge connections
    edge_attr = []   # List to store edge attributes (distance, angles)
    
    # Iterate over each particle in the DataFrame
    for _, row in df.iterrows():
        particle_id = row['parID']  # Unique particle identifier
        x0, y0, z0 = row['x0'], row['y0'], row['z0']  # Particle coordinates
        
        # Process N nearest neighbors
        for i in range(1, N + 1):
            neighbour_id = row[f'n{i}']  # Neighbor identifier
            x_neighbor, y_neighbor, z_neighbor = row[f'x{i}'], row[f'y{i}'], row[f'z{i}']
            
            # Add edge (particle -> neighbor)
            edge_index.append([particle_id, neighbour_id])
            
            # Compute relative position vector
            relative_position = torch.tensor([x_neighbor - x0, y_neighbor - y0, z_neighbor - z0], dtype=torch.float64)
            
            # Compute spherical coordinates
            r = relative_position.norm()  # Distance between particles
            phi = torch.atan2(relative_position[1], relative_position[0]).unsqueeze(0)  # Azimuthal angle
            theta = torch.acos(relative_position[2] / r).unsqueeze(0)  # Polar angle
            
            # Store edge attributes (distance, azimuthal angle, polar angle)
            edge_attr.append([r, phi, theta])
    
    # Convert lists to tensors
    edge_index = torch.tensor(edge_index, dtype=torch.long).t()  # Shape: [2, num_edges]
    edge_attr = torch.tensor(edge_attr, dtype=torch.float64)  # Shape: [num_edges, 3]
    
    # Normalize forces using the Reynolds number factor
    if force=='Cd0':
        forces = torch.tensor((df[[force]].values - df[force].mean()) * RE, dtype=torch.float64).to(device)
    else:
        forces = torch.tensor((df[[force]].values * RE), dtype=torch.float64).to(device)
        
    # Extract and store rotation angles as node features
    rotation_angles = torch.tensor(df2[['phi', 'psi']].values, dtype=torch.float64).to(device)
    
    # Construct PyTorch Geometric Data object
    data = Data(x=rotation_angles, edge_index=edge_index, edge_attr=edge_attr, y=forces)
    
    return data

def create_masks(data, train_ratio=0.8, seed=None):
    """
    Create training and validation masks for a graph dataset.

    Args:
        data (torch_geometric.data.Data): Graph data object.
        train_ratio (float, optional): Fraction of nodes used for training. Default is 0.8.
        seed (int, optional): Random seed for reproducibility. Default is None.

    Returns:
        torch_geometric.data.Data: Updated data object with `train_mask` and `val_mask` attributes.
    """
    num_nodes = data.num_nodes  # Total number of nodes in the graph
    indices = list(range(num_nodes))  # Create a list of node indices

    # Set random seed for reproducibility
    if seed is not None:
        random.seed(seed)
    
    # Shuffle the node indices
    random.shuffle(indices)

    # Split into training and validation sets
    train_size = int(train_ratio * num_nodes)
    train_indices = indices[:train_size]
    val_indices = indices[train_size:]

    # Initialize boolean masks
    data.train_mask = torch.zeros(num_nodes, dtype=torch.bool)
    data.val_mask = torch.zeros(num_nodes, dtype=torch.bool)

    # Assign mask values based on the split
    data.train_mask[train_indices] = True
    data.val_mask[val_indices] = True

    return data
