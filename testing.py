import mol2class
from tqdm import tqdm
import os
import torch
import torch.nn as nn
import torch.optim as optim
import sys


class GNN(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(GNN, self).__init__()
        self.fc1 = torch.nn.Linear(input_dim, hidden_dim)
        self.fc2 = torch.nn.Linear(hidden_dim, output_dim)

    def forward(self, x, adj):
        x = self.fc1(x)
        x = torch.relu(x)
        x = self.fc2(x)
        x = torch.sigmoid(x)
        return x


# Step 2: Initialize the GNN model
input_dim = 22
hidden_dim = 128
output_dim = 1
gnn_model = GNN(input_dim, hidden_dim, output_dim)

rootdir = 'scPDB/'
file = '4f0e_1'

# Load the saved model
gnn_model.load_state_dict(torch.load('gnn_model.pth'))

mol = mol2class.readMol2(rootdir+file+'/protein.mol2') # Get the molecule into the readMol2 class
numAtoms = mol.numAtoms()
for i in tqdm(range(numAtoms), desc="Generating Matrices...", file=sys.stdout):
    adj_matrix = mol.adjacencyMatrix(mol.atoms()[i])
    molSol = mol2class.Mol2ligand(rootdir+file+'/cavity6.mol2') # Get the molecule into the readMol2 class
    feature_matrix = molSol.SolutionsFeatureMatrix(mol.featureMatrix(mol.atoms()[i]))
    adj_matrix = torch.tensor(adj_matrix)
    feature_matrix = torch.tensor(feature_matrix)
    # Step 5: Train the GNN model
    adj_matrix = adj_matrix.to(dtype=torch.float32)
    feature_matrix = feature_matrix.to(dtype=torch.float32)
    output = gnn_model(feature_matrix, adj_matrix)

# For each position we can compute the probability of being a binding site by the mean of the outputs for a given position