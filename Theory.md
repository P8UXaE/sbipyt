
# Background and cientific explanation 

SASA (Solvent Accessible Surface Area) is a program that calculates the solvent-accessible surface area of a molecule, which is a measure of the exposed surface area of a molecule that is available for interaction with solvent molecules. The SASA program calculates the surface area of a molecule by using an algorithm that defines the solvent-accessible surface as the surface that can be reached by a probe sphere that rolls along the surface of the molecule without penetrating it.
The output of the SASA program typically includes a feature matrix and an adjacent matrix that describe the molecule's surface area and connectivity, respectively.
The feature matrix contains a row for each atom in the molecule and a column for each feature, such as the SASA, atomic radius, or charge. The SASA value for each atom is usually included as one of the features in the matrix.

The adjacent matrix describes the connectivity between atoms in the molecule, typically in the form of a sparse matrix where each row and column corresponds to an atom, and the entries indicate the strength of the bond between the atoms. In some cases, the adjacent matrix may also include information on non-covalent interactions, such as hydrogen bonds or Van der Waals interactions.