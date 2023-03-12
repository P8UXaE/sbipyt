import mol2class
from tqdm import tqdm
import os
import torch
import torch.nn as nn
import torch.optim as optim



# Step 1: Define the GNN model architecture
# class GNN(nn.Module):
#     def __init__(self, input_dim, hidden_dim, output_dim):
#         super(GNN, self).__init__()
#         self.fc1 = nn.Linear(input_dim, hidden_dim)
#         self.fc2 = nn.Linear(hidden_dim, output_dim)
#         self.act = nn.ReLU()

#     def forward(self, feature_matrix, adj_matrix):
#         x = torch.cat((feature_matrix, adj_matrix), dim=1)
#         x = self.fc1(x)
#         x = self.act(x)
#         x = self.fc2(x)
#         return x

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

# Step 3: Define the loss function and optimizer
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(gnn_model.parameters(), lr=0.001)



rootdir = 'scPDB/'

# Step 4: Loas the data for training
for file in os.listdir(rootdir): # This will get all the proteins inside the scPDB folder
    print('├─'+file)
    print('│ ├─'+'protein.mol2')
    mol = mol2class.readMol2(rootdir+file+'/protein.mol2') # Get the molecule into the readMol2 class
    numAtoms = mol.numAtoms()
    for i in tqdm(range(numAtoms), desc="Generating Matrices..."):
        adj_matrix = mol.adjacencyMatrix(mol.atoms()[i])
        molSol = mol2class.Mol2ligand(rootdir+file+'/cavity6.mol2') # Get the molecule into the readMol2 class
        feature_matrix, solutions = molSol.SolutionsFeatureMatrix(mol.featureMatrix(mol.atoms()[i]))

        adj_matrix = torch.tensor(adj_matrix)
        feature_matrix = torch.tensor(feature_matrix)
        
        # print(adj_matrix.dtype, feature_matrix.dtype)
        # print(adj_matrix.shape,feature_matrix.shape)
        # print(adj_matrix,feature_matrix)
        # x = torch.cat((feature_matrix.to(torch.float64), adj_matrix.to(torch.float64)), dim=1)
        # print(x)
        # Step 5: Train the GNN model
        adj_matrix = adj_matrix.to(dtype=torch.float32)
        feature_matrix = feature_matrix.to(dtype=torch.float32)
        
        output = gnn_model(feature_matrix, adj_matrix)
        
        solutions2 = []
        for i in solutions:
            solutions2.append([i])
        solutions = solutions2
        solutions = torch.tensor(solutions)
        solutions = solutions.to(dtype=torch.float32)

        # print(output, solutions)

        loss = criterion(output, solutions)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        # Step 6: Save the GNN model
        torch.save(gnn_model.state_dict(), 'gnn_model.pth')
        # Step 7: Load the GNN model
        gnn_model.load_state_dict(torch.load('gnn_model.pth'))

    print('│ ├─'+'cavity6.mol2')




# # Step 4: Load the data for training
# for batch_idx, (feature_matrix, adj_matrix, labels) in enumerate(train_loader):
#     # Step 5: Train the GNN model
#     output = gnn_model(feature_matrix, adj_matrix)
#     loss = criterion(output, labels)
#     optimizer.zero_grad()
#     loss.backward()
#     optimizer.step()

#     # Step 6: Save the GNN model
#     torch.save(gnn_model.state_dict(), 'gnn_model.pth')

#     # Step 7: Load the GNN model
#     gnn_model.load_state_dict(torch.load('gnn_model.pth'))




# def incrementalGNN(rootdir):
#     for file in os.listdir(rootdir): # This will get all the proteins inside the scPDB folder
#         print('├─'+file)

#         print('│ ├─'+'protein.mol2')
#         mol = mol2class.readMol2(rootdir+file+'/protein.mol2') # Get the molecule into the readMol2 class
#         numAtoms = mol.numAtoms()
#         for i in tqdm(range(numAtoms), desc="Generating Matrices..."):

#             adjacencyMatrix = mol.adjacencyMatrix(mol.atoms()[i])
#             molSol = mol2class.Mol2ligand(rootdir+file+'/cavity6.mol2') # Get the molecule into the readMol2 class
#             featureMatrix, solutions = molSol.SolutionsFeatureMatrix(mol.featureMatrix(mol.atoms()[i]))


#         print('│ ├─'+'cavity6.mol2')

# incrementalGNN('scPDB/')

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