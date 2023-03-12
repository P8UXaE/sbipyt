import mol2class
from tqdm import tqdm
import os
import numpy as np



def create_graphs(mol_file, feat_dim=16):
    """
    Creates adjacency matrices and feature matrices for each atom in a mol2 file, where each atom is represented as a graph
    with its 16 neighboring atoms.
    Args:
        mol2_file (str): Path to the mol2 file.
        feat_dim (int): Number of features to include for each atom. Default is 16.
    Returns:
        adjacency_matrices (list of np.ndarray): List of adjacency matrices for each atom.
        feature_matrices (list of np.ndarray): List of feature matrices for each atom.
    """

    mol = mol2class.readMol2(mol_file) # Initializate readMol2 class
    numAtoms = mol.numAtoms()
    adjacency_matrices = []
    feature_matrices = []
    for i in tqdm(range(numAtoms), desc="Generating Matrices..."):
        # print(mol.getNeighbors(mol.atoms()[i], k=16))
        adjacency_matrices.append(mol.adjacencyMatrix(mol.atoms()[i]))
        # print('-'*100)
        feature_matrices.append(mol.featureMatrix(mol.atoms()[i]))
        # for i in mol.featureMatrix(mol.atoms()[i]):
        #     print(i)
        # break
    # for i, j in zip(mol.atoms(), mol.sasa()):
    #     print(i, j)
    # for i in mol.sasa():
    #     print(i)




    return adjacency_matrices, feature_matrices




# print(create_graphs('scPDB/1a2b_1/protein.mol2'))

def incrementalGNN(rootdir):

    for file in os.listdir(rootdir): # This will get all the proteins inside the scPDB folder
        print('├─'+file)

        print('│ ├─'+'protein.mol2')
        mol = mol2class.readMol2(rootdir+file+'/protein.mol2') # Get the molecule into the readMol2 class
        numAtoms = mol.numAtoms()
        for i in tqdm(range(numAtoms), desc="Generating Matrices..."):
            # print(mol.featureMatrix(mol.atoms()[i]))
            adjacencyMatrix = mol.adjacencyMatrix(mol.atoms()[i])
            molSol = mol2class.Mol2ligand(rootdir+file+'/cavity6.mol2') # Get the molecule into the readMol2 class
            featureMatrix = molSol.SolutionsFeatureMatrix(mol.featureMatrix(mol.atoms()[i]))


            # print(molSol)

        print('│ ├─'+'cavity6.mol2')
        #molSol = mol2class.Mol2ligand(rootdir+file+'/cavity6.mol2', mol) # Get the molecule into the readMol2 class
                


incrementalGNN('scPDB/')