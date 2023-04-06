#!/usr/bin/python
import sys
import argparse
import torch
import gnnclass as gnc
import proteinclass
from tqdm import tqdm
import numpy as np





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
    input_dim = 57
    hidden_dim = 4
    output_dim = 1
    gnn_model = gnc.GNN(input_dim, hidden_dim, output_dim)
    gnn_model.load_state_dict(torch.load('gnn_model.pth'))


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
        adj_matrix = mol.adjacencyMatrix(mol.atoms()[i])
        feature_matrix = mol.featureMatrix(mol.atoms()[i])
        adj_matrix = torch.tensor(adj_matrix)
        feature_matrix = torch.tensor(feature_matrix)
        adj_matrix = adj_matrix.to(dtype=torch.float32)
        feature_matrix = feature_matrix.to(dtype=torch.float32)
        # print(adj_matrix, feature_matrix)
        output = gnn_model(feature_matrix, adj_matrix)
        c = 0
        for calc, neigh, sol in zip(feature_matrix.tolist(), mol.getNeighbors(mol.atoms()[i]), output.tolist()):
            neigh = neigh[0]
            sas = calc[1]
            sol = sol[0]
            solutions_probability.setdefault(neigh[6], {})
            solutions_probability[neigh[6]].setdefault(neigh[0], {})
            solutions_probability[neigh[6]][neigh[0]]['sasa'] = sas
            solutions_probability[neigh[6]][neigh[0]]['solution'+str(c)] = sol
            solutions_probability[neigh[6]][neigh[0]]['position'] = [neigh[2], neigh[3], neigh[4]]
            c += 1
            
            
    print(solutions_probability)
    