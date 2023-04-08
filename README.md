# SBI-PYTHON project
## BITER: a Binding sITEs detectoR for proteins

The main approach to detect binding sites on proteins is by ML. It also has a geometric approach that is able to detect pockets. The optimal version of python to run the algorithm is v3.9 because is the most suitable to work with packages such as pythorch.


## Tutorial

It is important to have completed those steps:
1.	Have downloaded the biter folder
2.	Have executed the setup.py in order to create the necessary environment

To run biter you just have to type the input file type and the file itself (pdb or mol2). It will generate a solution file that can be called using chimera to show the distribution of probabilities through all the residues to take part on a binding site.
To see all the options you can type and mer information
```
$ python biter.py --help
```

By typing
```
$ python biter.py -i pdb 1mee.pdb 
```
The program will run all, the ML approach and the geometric approach. This will generate a pdb file that is a single-chain pdb, more suitable to show the results in chimera. As will be printed in the terminal you can run 3 different commands:
```
$ chimera 1mee.pdb 1mee_pocketPoints.pdb
$ chimera 1mmee_chimera.cmd
$ chimera 1mmee_chimera.cmd 1mee_pocketPoints.pdb
```
The first command will open a chimera window showing the points found that are inside a protein pocket. Then, you can go to Favorites > Command Line to activate chimera's command line and type:
```
sel: X
```
Where X is a residue number. In the console you will see which are the top 3 clusters of those pocket points. Also, all clusters ordered by distance between the points will appear.
The second command will pop up a chimera window showing the probability of each residue of taking part in a binding process. The theory of everyting is lately explained. The darker the color, the higher the probability.
```
PERCENTAGE  RGB 0-1
0-10%       1,0.96,0.9
10-20%      1,0.91,0.8
20-30%      1,0.85,0.65
30-40%      1,0.75,0.47
40-50%      1,0.66,0.3
50-60%      1,0.57,0.17
60-70%      0.99,0.49,0.08
70-80%      0.97,0.4,0.03
80-90%      0.91,0.35,0.05
90-100%     0.85,0.28,0.06
```
Finnaly, the third command will open a chimera window showing both, the ML approach and the geometric one.  
Here some possible commands are shown:  
```
$ python biter.py -i pdb -b protein.pdb myprotein1
```
This calls the program using a pdb type file. It will only generate the ML solution (-b stands for biter, so you can also type --biter). The output files will have 'myprotein1' name on it instead of 'protein'.
```
$ python biter.py -i mol2 -p 1iki.mol2
```
This calls the program using a mol2 type file. It will only generate the geometric solution (-p stands for pocket, so you can also type --pocket). The output files will have '1iki' name on it as a second argument has not been passed to the command.

You can find some examples of the results and the commands used in the examples folder.


## Training Set
The training set is obtained from *http://bioinfo-pharma.u-strasbg.fr/scPDB/*. It is a mannualy curated database that contains proteins with known binding sites.

This tree represents how files are stored inside the scPDB folder

```
scPDB
├── 1a2b_1
│   ├── IFP.txt
│   ├── cavity6.mol2
│   ├── cavityALL.mol2
│   ├── ints_M.mol2
│   ├── ligand.mol2
│   ├── ligand.sdf
│   ├── protein.mol2
│   └── site.mol2
├── 1a2n_1
│   ├── IFP.txt
│   ├── cavity6.mol2
│   ├── cavityALL.mol2
│   ├── ints_M.mol2
│   ├── ligand.mol2
│   ├── ligand.sdf
│   ├── protein.mol2
│   └── site.mol2
└── 1a4r_1
    ├── IFP.txt
    ├── cavity6.mol2
    ├── cavityALL.mol2
    ├── ints_M.mol2
    ├── ligand.mol2
    ├── ligand.sdf
    ├── protein.mol2
    └── site.mol2
```

The protein structure is extracted from the protein.mol2 file. The atoms that are binding atoms (ommiting H) are labeled thanks to both cavity and site. Cavity is just a mol2 file with spaced dots that represent the cavity where the ligand is placed. The site is a residue subset that are near the cavity. Using those and a discance of x Å we can determine which atoms should be labeled as bindin atoms inside the cavity. Not taking into account any type of bond or Van der Waals force.
The training process:
1. Open the .pth file (if exists) that contains the model.
2. Using the model, do the calculus and the parameter modification of the model for every atom in a protein.
3. Once all the protein has been evaluated and has been used in order to improve the model, the protein folder is written down in a list in order to not be used again to train the model, given that the program is coded in order to use every protein in scPDB folder to train the model.

## Create working environment

If the setup.py file does not work, you can do it manually. First of all you need to create a virtual environment. In order to do this you can run the following comands in the parent folder.
```bash
$ python3.9 -m venv python3_9venv
```
You can activate the virtual environment by typing the following
```bash
$ source python3_9venv/bin/activate
```
And you can deactivate the environment by
```bash
$ deactivate
```
To remove the virtual env:
```bash
$ rm -rf python3_9venv
```
Once you are inside the environment, you can use pip to install any package
```bash
$ pip install numpy
```


Check modules:
$ pip list
$ pip3 list
$ pip freeze > requirements.txt

Remove module:
$ pip remove <package-name>

# setup 
Used to package all the modules or application for the project into a distributable package that can be installed by other users. The script contains information about the package, including its name, version, author, description, and dependencies, among other things.
the setuptools package is needed 
$ pip install setuptools

load all the requirements for the project 
$ python Setup.py install
$ python Setup.py build

model training 


## PDB style


The pdb file must be like the following example in order to be integrated in the python class:
(1mee - *https://www.rcsb.org/structure/1MEE*)
```bash
ATOM      1  N   ALA A   1     -16.582  25.909  46.648  1.00 25.87           N  
ATOM      2  CA  ALA A   1     -16.120  24.544  46.323  1.00 60.37           C  
ATOM      3  C   ALA A   1     -15.090  24.573  45.193  1.00 39.03           C  
ATOM      4  O   ALA A   1     -15.261  25.292  44.195  1.00 25.75           O  
ATOM      5  CB  ALA A   1     -17.257  23.611  45.965  1.00 31.78           C  
ATOM      6  N   GLN A   2     -14.060  23.750  45.366  1.00 24.97           N  
ATOM      7  CA  GLN A   2     -13.012  23.760  44.315  1.00 17.40           C  
ATOM      8  C   GLN A   2     -13.100  22.564  43.399  1.00 15.42           C  
ATOM      9  O   GLN A   2     -13.394  21.452  43.865  1.00 14.10           O  
ATOM     10  CB  GLN A   2     -11.668  23.893  45.007  1.00 14.39           C  
ATOM     11  CG  GLN A   2     -10.450  23.827  44.068  1.00 12.20           C  
```
This pdb type won't work:
(8ad1 - *https://www.rcsb.org/structure/8AD1*)
```bash
ATOM    433  P    DC N 263     157.725 135.829 113.933  1.00255.48           P  
ATOM    434  OP1  DC N 263     158.817 136.674 113.402  1.00255.48           O  
ATOM    435  OP2  DC N 263     156.489 135.643 113.142  1.00255.48           O  
ATOM    436  O5'  DC N 263     158.316 134.384 114.274  1.00255.48           O  
ATOM    437  C5'  DC N 263     159.510 134.274 115.034  1.00255.48           C  
ATOM    438  C4'  DC N 263     159.692 132.864 115.565  1.00255.48           C  
ATOM    439  O4'  DC N 263     158.576 132.513 116.426  1.00255.48           O  
ATOM    440  C3'  DC N 263     159.760 131.766 114.501  1.00255.48           C  
ATOM    441  O3'  DC N 263     160.725 130.801 114.890  1.00255.48           O  
ATOM    442  C2'  DC N 263     158.349 131.183 114.538  1.00255.48           C  
ATOM    443  C1'  DC N 263     158.064 131.264 116.024  1.00255.48           C  
```

# 