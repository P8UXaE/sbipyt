#!/usr/bin/python
import sys
import argparse
import torch
import gnnclass as gnc
import proteinclass
from tqdm import tqdm
import numpy as np
import torch.nn as nn





if __name__=='__main__':
    ######################################################
            ## THIS IS THE COMMAND LINE OPTIONS ##
    ######################################################
    parser = argparse.ArgumentParser(description="Binding sITEs detectoR (BITER) is a GNN model that predicts the probability of the different regions of a protein to be binding site. Example: $ python biter.py -i mol2 protein_example_biter.mol2")
    parser.add_argument('-i', '--input',
                        dest= "intype",
                        action= None,
                        default= "pdb",
                        choices= ["pdb", "mol2"],
                        help= "Input file type [pdb, mol2], default pdb")
    
    parser.add_argument('fileRoute',
                        help= "Route of the file to be analyzed")

    options = parser.parse_args()
    # print(options.intype)
    # print(options.fileRoute)


    ######################################################
            ## INITIALIZATE GNN MODEL ##
    ######################################################
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
    
    mlp_model.load_state_dict(torch.load('gnn_model2.pth'))


    ######################################################
            ## THIS IS THE MOL2 OPTION ##
    ######################################################
    if options.intype == 'mol2':
        mol = proteinclass.readMod2(options.fileRoute) # Get the molecule into the readMol2 class

    ######################################################
            ## THIS IS THE PDB OPTION ##
    ######################################################
    if options.intype == 'pdb':
        mol = proteinclass.readpdb(options.fileRoute)

    

    ######################################################
            ## GET THE SOLUTION FOREACH ATOM ##
    ######################################################
    solutions_probability = {} # Structure: residue - atom - sasa, solution
    numAtoms = mol.numAtoms()
    for i in tqdm(range(numAtoms), desc="Generating Matrices...", file=sys.stdout):
        feature_matrix = mol.featureMatrix(mol.atoms()[i])
        feature_matrix = [item for feature in feature_matrix for item in feature]
        feature_matrix = torch.tensor(feature_matrix, requires_grad=True)
        output = mlp_model(feature_matrix)
        output = output.item()

        # print(mol.getNeighbors(mol.atoms()[i]))

        residue = mol.getNeighbors(mol.atoms()[i])[0][0][6]
        atom = mol.getNeighbors(mol.atoms()[i])[0][0][0]
        x = mol.getNeighbors(mol.atoms()[i])[0][0][2]
        y = mol.getNeighbors(mol.atoms()[i])[0][0][3]
        z = mol.getNeighbors(mol.atoms()[i])[0][0][4]

        solutions_probability.setdefault(residue, {})
        solutions_probability[residue].setdefault(atom, {})

        solutions_probability[residue][atom]['sasa'] = feature_matrix.tolist()[1]
        solutions_probability[residue][atom]['solution'] = output
        solutions_probability[residue][atom]['position'] = [x, y, z]
            
    print(solutions_probability)
    