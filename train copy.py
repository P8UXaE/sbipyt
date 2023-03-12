import mol2class
from tqdm import tqdm
import os
import torch
import torch.optim as optim


def incrementalGNN(rootdir):
    for file in os.listdir(rootdir): # This will get all the proteins inside the scPDB folder
        print('├─'+file)

        print('│ ├─'+'protein.mol2')
        mol = mol2class.readMol2(rootdir+file+'/protein.mol2') # Get the molecule into the readMol2 class
        numAtoms = mol.numAtoms()
        for i in tqdm(range(numAtoms), desc="Generating Matrices..."):

            adjacencyMatrix = mol.adjacencyMatrix(mol.atoms()[i])
            molSol = mol2class.Mol2ligand(rootdir+file+'/cavity6.mol2') # Get the molecule into the readMol2 class
            featureMatrix, solutions = molSol.SolutionsFeatureMatrix(mol.featureMatrix(mol.atoms()[i]))


        print('│ ├─'+'cavity6.mol2')

incrementalGNN('scPDB/')

            # if 1 in solutions:
            #     print('There is a solution')
            #     print(solutions)

'''
import os
import torch
import torch.optim as optim
import mol2class

# Define the GNN model
class GNN(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(GNN, self).__init__()
        self.conv1 = torch.nn.Linear(input_dim, hidden_dim)
        self.conv2 = torch.nn.Linear(hidden_dim, output_dim)

    def forward(self, adj, feat):
        h = torch.relu(self.conv1(torch.matmul(adj, feat)))
        out = self.conv2(torch.matmul(adj, h))
        return out

# Initialize the GNN model
input_dim = 64 # dimension of feature matrix
hidden_dim = 128 # dimension of hidden layer
output_dim = 1 # dimension of output layer
gnn_model = GNN(input_dim, hidden_dim, output_dim)

# Define the optimizer
learning_rate = 0.001
optimizer = optim.Adam(gnn_model.parameters(), lr=learning_rate)

# Train the model incrementally
rootdir = "scPDB/"

for file in os.listdir(rootdir):
    print('├─'+file)
    print('│ ├─'+'protein.mol2')
    
    # Load the protein molecule
    mol = mol2class.readMol2(rootdir+file+'/protein.mol2')
    numAtoms = mol.numAtoms()
    
    # Load the ligand molecule
    molSol = mol2class.Mol2ligand(rootdir+file+'/cavity6.mol2')
    
    for i in range(numAtoms):
        # Generate adjacency matrix
        adj = mol.adjacencyMatrix(mol.atoms()[i])
        adj_tensor = torch.from_numpy(adj).float()

        # Generate feature matrix
        feat = molSol.SolutionsFeatureMatrix(mol.featureMatrix(mol.atoms()[i]))
        feat_tensor = torch.from_numpy(feat).float()

        # Compute predicted output
        output = gnn_model(adj_tensor, feat_tensor)
        
        # Compute ground truth output (replace with your own implementation)
        ground_truth = 0.5
        
        # Compute loss and update parameters
        loss = torch.nn.MSELoss()(output, ground_truth)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
    print('│ ├─'+'cavity6.mol2')
    
# Evaluate the trained model (replace with your own implementation)
val_loss = 0.1
test_loss = 0.2
'''