import os

rootdir = 'scPDB/'
for file in os.listdir(rootdir): # This will get all the proteins inside the scPDB folder
    print('├─'+file)
    for file2 in os.listdir(rootdir+str(file)): # This will get all the files inside each protein file
        if file2 == 'protein.mol2' or file2 == 'ligand.mol2':
            print('│ ├─'+file2)
