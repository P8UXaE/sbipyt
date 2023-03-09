# sbipyt
SBI-PYT project

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

