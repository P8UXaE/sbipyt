#!/usr/bin/python
import sys
import argparse
import torch
import proteinclass
from tqdm import tqdm
import numpy as np
import torch.nn as nn
import os
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import pdist, squareform




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
                        help= "input file type [pdb, mol2], default pdb")
    
    parser.add_argument('-b', '--biter',
                        action= 'store_true',
                        help= "just run the ML (geometric integred) program")
    
    parser.add_argument('-p', '--pocket',
                        action= 'store_true',
                        help= "just run the geometric pocket point finder")
    
    parser.add_argument('fileRoute',
                        help= "route of the file to be analyzed")

    options = parser.parse_args()
    # print(options.intype)
    protein_name = os.path.splitext(os.path.basename(options.fileRoute))[0]


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
            ## RUN THE GEOMETRIC PROGRAM ##
    ######################################################
    if options.pocket and not options.biter or not options.pocket and not options.biter:
        bs = open(protein_name+'_pocketPoints.pdb', 'w')
        n = 1
        pocket_points = mol.pocketList()

        # for i in pocket_points:
        #     bs.write("ATOM"+" "*(7-len(str(n)))+str(n)+"  O   HOH X"+" "*(4-len(str(n)))+str(n)+" "*(8-len(str(str(i[0]).split('.')[0])))+str(format(i[0], ".3f"))+" "*(4-len(str(str(i[1]).split('.')[0])))+str(format(i[1], ".3f"))+" "*(4-len(str(str(i[2]).split('.')[0])))+str(format(i[2], ".3f"))+"  1.00 30.00           O  \n")
        #     n += 1

        ### CLUSTERING ###
        # Define a range of cluster numbers to test
        cluster_range = range(2, 30)
        # Calculate the silhouette score for each number of clusters
        silhouette_scores = []
        for n_clusters in cluster_range:
            kmeans = KMeans(n_clusters=n_clusters, n_init=3, init='k-means++', random_state=42)
            kmeans.fit(pocket_points)
            labels = kmeans.labels_
            silhouette_scores.append(silhouette_score(pocket_points, labels))
        maximums = []
        optimal_num_clusters = 0
        # Extract the optimal number of clusters
        for i in range(len(silhouette_scores)-1):
            if silhouette_scores[i+1]<silhouette_scores[i]:
                maximums.append(i)
                if len(maximums) > 3:
                    optimal_num_clusters = i+2
                    break
        if optimal_num_clusters == 0:
            optimal_num_clusters = maximums[-1]+2
        
        # print(optimal_num_clusters)

        # Fit the KMeans model to the data
        kmeans = KMeans(n_clusters=optimal_num_clusters, n_init=3, init='k-means++', random_state=42)
        kmeans.fit(pocket_points)
        # Get the cluster labels for each point
        labels = kmeans.labels_
        # Create a list of points that belong to each cluster
        clusters = [[] for _ in range(optimal_num_clusters)]
        for i, label in enumerate(labels):
            clusters[label].append(pocket_points[i])
        pocket_atoms = ['OG', 'NZ', 'CZ', 'HB']
        pocket_residues = ['SER', 'LYS', 'PHE', 'ASP']
        # Print the points that belong to each cluster
        for j, cluster in enumerate(clusters):
            if j > 3:
                color = j%4
            else:
                color = j
            for i in cluster:
                bs.write("ATOM"+" "*(7-len(str(n)))+str(n)+" "+pocket_atoms[color]+"   "+pocket_residues[color]+" X"+" "*(4-len(str(j)))+str(j)+" "*(8-len(str(str(i[0]).split('.')[0])))+str(format(i[0], ".3f"))+" "*(4-len(str(str(i[1]).split('.')[0])))+str(format(i[1], ".3f"))+" "*(4-len(str(str(i[2]).split('.')[0])))+str(format(i[2], ".3f"))+"  0.00 00.00           "+pocket_atoms[color][0]+"  \n")
                n += 1
            # print(f"Cluster {i}: {cluster}")
        # Calculate the average distance between points within each cluster
        avg_distances = []
        for cluster in clusters:
            distances = squareform(pdist(cluster))
            avg_distance = distances.sum() / (len(cluster) * (len(cluster) - 1))
            avg_distances.append(avg_distance)
        # Find the cluster with the lowest average distance
        best_cluster_index = avg_distances.index(min(avg_distances))
        # best_cluster = clusters[best_cluster_index]
        print('-'*30)
        print('Please, check the top 3 clusters')
        print('Go to Chimera, open the protein requested and open the _pocketPoints.pdb file. You can use:')
        print('$ chimera '+options.fileRoute+" "+protein_name+'_pocketPoints.pdb')
        print('Then activate the command line; Favorites > Command Line')
        print('Type ´sel: X´ to select and visualize the Clusters detected of pocket binding sites.')
        print('-'*30)
        print('Best cluster distance:', best_cluster_index)
        top3_indices = sorted(range(len(avg_distances)), key=lambda i: avg_distances[i])[:3]
        print('Top 3 clusters distance:', top3_indices)
        print('All clusters:', sorted(range(len(avg_distances)), key=lambda i: avg_distances[i]))
        print('-'*30)

        # # Plot the silhouette scores
        # plt.plot(cluster_range, silhouette_scores)
        # plt.xlabel('Number of Clusters')
        # plt.ylabel('Silhouette Score')
        # plt.title('Silhouette Analysis')
        # plt.show()

        bs.write("END                                                                             \n")
    
    ######################################################
            ## GET THE SOLUTION FOREACH ATOM ##
    ######################################################
    if not options.pocket and not options.biter or options.biter and not options.pocket:
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
    