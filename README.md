# SBI-PYTHON project
## BITER: a Binding sITEs detectoR for proteins

*https://github.com/P8UXaE/sbipyt*

Authors: Pau Pujol Vives & Junghua Ye

Contact: paupujolvives@gmail.com junhuay00@gmail.com

## Introduction

The detection of binding sites is an important field of study. The development of an algorithm able to detect possible binding sites in several and different molecules is necessary in order to have more specific and personalized drugs and to keep increassing the knowledge level in protein-protein interactions and protein-ligand interactions.
Biter integrates both, geometric pocket detection and machine learning (ML) algorithms in order to get closer to the mentioned objective.

## Tutorial

It is important to have completed those steps:
1.	Have downloaded the biter folder
2.	Have executed the setup.sh in order to create the necessary environment

Before running the biter.py script you must have installed all the necessary packages, either in the environment or in the main pyhton libraries folder. If you are using the environment you have to activate it.
```
$ source python3_9venv/bin/activate
```
This command will turn on the environment where you will have executed the setup.py file, and consequently installed the necessary packages. Once you have used the tool, you can deactivate the environment.
```
$ deactivate
```

To run biter you have to type the input file type and the file itself (pdb or mol2). It will generate a solution file that should be called using chimera to show the distribution of probabilities through all the residues that take part on a binding site and the pocket detected with the geometric algorithm.
To see all the options you can type and help and get more information (run inside biter/)
```
$ python biter.py --help
```
### Running the program
```
$ python biter.py -i pdb 1mee.pdb 
```
The program will run all, the ML approach and the geometric approach. This will generate a pdb file that is a single-chain pdb (only happens when working with pdb file type), more suitable to show the results in chimera. As will be printed in the terminal you can run 3 different commands:
```
$ chimera 1mee.pdb 1mee_pocketPoints.pdb
$ chimera 1mmee_chimera.cmd
$ chimera 1mmee_chimera.cmd 1mee_pocketPoints.pdb
```
The first command will open a chimera window showing the points found that are inside a protein pocket. Then, you can go to Favorites > Command Line to activate chimera's command line. Then type:
```
sel: X
```
Where X is a cluster number. In the terminal where you have executed the biter program will appear which are the top 3 clusters of those pocket points. Also, all clusters ordered by distance between the points will appear.
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
### Comands examples
Here some possible commands are shown:  
```
$ python biter.py -i pdb -b protein.pdb myprotein1
```
This calls the program using a pdb type file. It will only generate the ML solution (-b stands for biter, so you can also type --biter). The output files will have 'myprotein1' name instead of 'protein'.
```
$ python biter.py -i mol2 -p 1iki.mol2
```
This calls the program using a mol2 type file. It will only generate the geometric solution (-p stands for pocket, so you can also type --pocket). The output files will have '1iki' name on it as the second argument has not been passed to the command.

You can find some examples of the results and the commands used in the examples folder.


## Automatic Setup 

Used to package all the modules or application for the project into a distributable package that can be installed by other users. The script contains information about the package, including its name, version, author, description, and dependencies, among other things. 

The setuptools package is needed.
```
$ pip install setuptools
```
You can run the setup.sh script to create a python working environment with a particular version of the python (3.9) and also install and load all the requirements for the project, you just have to activate the environment and you will be able to run the biter programm. 
Before make executable the setup.sh.

```
$ chmod +x setup.sh 
$ ./setup.sh
```


## Create working environment - Manual Setup

If the setup.sh file does not work or you prefer doing it yourself, you can do it manually. First of all you need to create a virtual environment. In order to do this you can run the following comands in the parent folder.
```
$ python3.9 -m venv python3_9venv
```
You can activate the virtual environment by typing the following
```
$ source python3_9venv/bin/activate
```
And you can deactivate the environment by
```
$ deactivate
```
To remove the virtual env (**CAREFUL**: only if you want to delete it)
```
$ rm -rf python3_9venv
```
Once you are inside the environment you can use pip to install any package, for example numpy
```
$ pip install numpy
```
Or run
```
$ python setup.py install
$ python setup.py build
```
To install all the packages.

## Training Set

The training set is obtained from *http://bioinfo-pharma.u-strasbg.fr/scPDB/*. It is a mannualy curated database that contains proteins (in a mol2 format) with known binding sites described with points.

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

The protein structure is extracted from the *protein.mol2* file. The atoms that belong to the protein (ommiting H) are labeled thanks to *cavity* file. It is a mol2 file with spaced dots that represent the cavity where the ligand is placed. Using this file and a discance of 3.5 Å we can determine which atoms should be labeled as binding atoms inside the cavity.

The pipeline training process:
1. Open the .pth file (if exists) that contains the model.
2. Go through atom by atom of the protein, predict whether it is a binding atom or not, and optimize the model.
3. Once all the protein has been evaluated and has been used in order to improve the model, the protein folder is written down in a list to not be used again to train the model, given that the program is coded to use every protein in scPDB folder to train the model.

This training method does not need any RAM memory to keep data from the training set, once the array of an atom has been used, it is replaced by the next one. In conventional techniques where you train a model with a dataframe, you need to keep big files with all the data, and the larger the data, the higher the RAM requirement.

In order to run the training algorithm you need to move to the training folder, and have downloaded the full scPDB data. Then just run
```
$ python train2.py 2> acc2.log
```
This will run the training script. You can see the accuracy only when the ML process has started, after the first 'Computing pocket points...'
```
$ tail -n 5 acc2.log
```
After training the model with 279 different molecules (found in already_trained_molecules.txt file) the accuracy is ~80%.

## PDB file working type

The pdb file must be like the following example in order to be integrated in the python class created:
(1mee - *https://www.rcsb.org/structure/1MEE*)
```
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
This pdb file won't work:
(8ad1 - *https://www.rcsb.org/structure/8AD1*)
```
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

## Theory

The binding sites detector is coded mainly as a machine learning approach. It first computes some data and using the pytorch package it trains an artificial intelligence that is capable to predict if each atom of a molecule takes part on a binding site or not. Lately, to filter the solutions, a clustering algorithm is computed in order to detect regions with high positive solution density. It wouldn't consider the Hidrogen molecules and train a model based on atomical level. To get every data row passed to the model, it must compute some features:
1. Every row contains the atom that we are evaluating and its 15 nearest neighbors.
2. For every neighbor:

    a) SASA value.

    b) Direction respect the main atom.

    c) The secondary structure of the residue of which the atom pertains.

    d) Lenard-Jones potential.

    e) Hydrophobicity of the residue of which the atom pertains.

    f) 10 nearest distances to geometric possible binding points.


### SASA value

The SASA (solvent-accessible surface area) value is computed using the Biopython package. It calculates the solvent-accessible surface area of a molecule, which is a measure of the exposed surface area of a molecule that is available for interaction with solvent molecules. The SASA program calculates the surface area by using an algorithm that defines the solvent-accessible surface as the surface that can be reached by a probe sphere that rolls along the surface of the molecule without penetrating it.
The output of the SASA program typically includes a feature matrix and an adjacent matrix that describe the molecule's surface area and connectivity, respectively.
The feature matrix contains a row for each atom in the molecule and a column for each feature, such as the SASA, atomic radius, or charge. The SASA value for each atom is usually included as one of the features in the matrix. [https://biopython.org/docs/dev/api/Bio.PDB.SASA.html]

### Direction respect the main atom

The objective of this precalculus is to both describe the protein structure using several features and adding as much relevant information as possible.
One of the ways to describe the structure is to know the direction of the 15 nearest neighbors of each atom.
This is basically computed using numpy arrays, and subtractions that result on a direction.

$$(AB)=A-B$$

### Secondary structure

The secondary structure of the protein is computed using the Φ and Ψ angles for the helix structure. For the beta sheet, the angle computed is the resulting between the N-H bond and the O.

Alpha helix:
For the alpha helix, the number of residues between one and another residue goes from 3 to 5. Also, the distance between those residues is less than 3.4 Å. Also Φ has to be 92 ± 35 %, and Ψ must be 98 ± 35 %.

$n_{rj} = n_{ri} + 3,4,5$

$d(n_{rj}, n_{ri}) < 3.4$

$92 - 32.2 < \Phi < 92 + 32.2$

$98 - 34.3 < \Psi < 98 + 34.3$


Beta sheet:
For the beta sheet first the distance between $N_i$ and $O_j$ must be < 3.2 Å. Then a H is placed in the same plane as $C_α-N-C$ at a distance of the midpoint of $C_α-C$, but in the opposite direction. Having placed the H in between $N_i$ and $O_j$ it is time to calculate the angle being H the vertex. The angle must be 180 ± 15 %.

### Lennard-Jones potential

It returns a matrix with the potential for every pair of atoms (of the 15 neighbors).

$$W(r)=4ε((σ/r)^{12}-(σ/r)^6 )$$

ε and σ are obtained from a library that has those values from each pair of atoms [https://github.com/choderalab/ambermini/blob/master/share/amber/dat/leap/parm/parm99.dat].

### Hydrophobicity

In order to transform the letter symbol of each residue and to add hydrophobicity information, a library that contains different hydrophobicity types is used [https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html#anote].

### Geometry based

To add more information about the geometry of the molecule we coded a script that calculates points around the protein that are placed on concave regions (pockets).

This script uses the SASA values to know which atoms do have accessible area, then it places points from [x-3, y-3, z-3] to [x+3, y+3, z+3] with a 3 value step around those atoms. Once it has the points it calculates whether each point is inside another atom or not. If not, it throws other points until reached a distance of 15 Å or a collision with another atom. If there is a collision it keeps whit the collision angles' values.

A draw in a 2D molecule representation is shown for a better understanding.

![Geometry approach 2d](geometryapproach.png "2D geometry approach algorithm simplifiaction visualisation")

In this 2D example we can see 2 atoms (red spheres) with SASA > 0. In this particular example each atom throws possible points (green dots). Those points can be inside of another atom. That would be the case of the green dots inside the protein. Later on, using the green points outside the protein, lines are thrown in 45º difference in all directions from distance 2 to 15 or until collision.

Once we have all the collisions we know the θ and ϕ angles ( in the case of 2D protein we would have just one angle) the surface of the source green point is calculated. If it is 0.5 or higher of the sphere (area equal or grater than 2π) we know we are in a concave site. That would be the case of 2 and 3 green dots.

In order to calculate the source sphere's area, it is divided into 4 tropics (that have a θ angle difference of 45º) and 8 meridians (that have a ϕ angle difference of 45º). The area of top and bottom poles is calculated with 3 existing points whether the 2 middle portions of the sphere (45º $\le$ θ $\le$ 135º).

The points that meet all this characteristics are then clustered and the higher density point regions are the output of the geometric algorithm.

Once the program has all the points, it adds 10 features foreach atom, those features being the 10 nearest distances to those points using KDTrees algorithm. 

## Further progress

While developing the algorithm some further investigation was thought:
1. Instead of working using single atoms, precalcualte the residue properties as a mean of all the atoms and use the residues to feed the ML model. This way less data to feed the ML is used and less variaty is seen.
2. Use techniques such as upgrading the minority group. This way we could rise the atoms or residues taking part in binding processes to higher levels, up to 50%. This could result in a better output and accuracy of the model.