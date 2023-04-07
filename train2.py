import proteinclass
from tqdm import tqdm
import os
import torch
import torch.nn as nn
import torch.optim as optim
import sys


class MLP(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, output_dim)
        self.act = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.fc1(x)
        x = self.act(x)
        x = self.fc2(x)
        x = self.sigmoid(x)
        return x


input_dim = 912
hidden_dim = 4
output_dim = 1
mlp_model = MLP(input_dim, hidden_dim, output_dim)

criterion = nn.BCELoss()
optimizer = optim.Adam(mlp_model.parameters(), lr=0.001)

rootdir = 'scPDB/'
accuracy = 0
iters = 1

with open('already_trained_molecules2.txt', 'r') as at:
    molecules_trained = at.read().splitlines()
at.close()

# Load the saved model
if os.path.isfile('gnn_model2.pth'):
    mlp_model.load_state_dict(torch.load('gnn_model2.pth'))

# Step 4: Loas the data for training
for file in os.listdir(rootdir): # This will get all the proteins inside the scPDB folder
    if file in molecules_trained:
        continue
    print('├─'+file)
    print('│ ├─'+'protein.mol2')
    mol = proteinclass.readMod2(rootdir+file+'/protein.mol2') # Get the molecule into the readMol2 class
    numAtoms = mol.numAtoms()
    for i in tqdm(range(numAtoms), desc="Generating Matrices...", file=sys.stdout):
        molSol = proteinclass.Mol2ligand(rootdir+file+'/cavity6.mol2') # Get the molecule into the readMol2 class

        solutions = molSol.SolutionsFeatureMatrix(mol.getNeighbors(mol.atoms()[i]), mol.sasaList())
        solution = solutions[0]
        # print('Solution:', solution)

        feature_matrix = mol.featureMatrix(mol.atoms()[i])
        feature_matrix = [item for feature in feature_matrix for item in feature]

        # print(feature_matrix)

        feature_matrix = torch.tensor(feature_matrix, requires_grad=True)
        # feature_matrix = feature_matrix.to(dtype=torch.float32)
        output = mlp_model(feature_matrix)
        output = output.item()
        # print('Output solution:', output)
        # solution = torch.tensor([solution])
        # solution = solution.to(dtype=torch.float32)
        # solution = torch.LongTensor(solution)

        # print(solution, output)


        accuracy += 1-abs(solution-output)
        print(str(accuracy/iters)+'\n',solution, output, file=sys.stderr)
        iters += 1
        
        # output = torch.tensor(output)
        # solution = torch.tensor(solution)


        loss = criterion(torch.tensor([output], dtype=torch.float32, requires_grad=True), torch.tensor([solution], dtype=torch.float32))
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

    # Step 6: Save the GNN model
    torch.save(mlp_model.state_dict(), 'gnn_model2.pth')
    # Step 7: Load the GNN model
    mlp_model.load_state_dict(torch.load('gnn_model2.pth'))
    print('│ ├─'+'cavity6.mol2')
    with open('already_trained_molecules2.txt', 'a') as at:
        at.write(file+'\n')
    at.close()
    molecules_trained.append(file)