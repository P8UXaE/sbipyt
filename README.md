# sbipyt
SBI-PYT project

Runs in python 3.9

There is a subset of 5 mol2 in scPDB folder. The .gz is 4..GB

When you run python3 gforatom.py it puts every molecule into the readMol2 class, where several features are extracted. Also it returns the sasa values, the adjacency matrixes...

## Training Set
The training set is obtained from *http://bioinfo-pharma.u-strasbg.fr/scPDB/*

This tree represents how the files are stored inside the scPDB folder

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

# git add . && git reset ./python3_9venv

# setup 
Used to package all the modules or application for the project into a distributable package that can be installed by other users. The script contains information about the package, including its name, version, author, description, and dependencies, among other things.
the setuptools package is needed 
$ pip install setuptools

load all the requirements for the project 
$ python Setup.py install

Create a virtual env:
$ python3.9 -m venv myenv (in order to install pytorch)
$ source myenv/bin/activate
$ deactivate
To remove the virtual env:
$ rm -rf myenv

Check modules:
$ pip list
$ pip3 list
$ pip freeze > requirements.txt

Remove module:
$ pip remove <package-name>
