import math
import numpy as np
import collections
from ctypes import *
from kdtrees import KDTree
import pyscf
# import sys
# sys.path.append('kdtrees.so')
# from kdtrees import KDTree

class readMol2():
    # __slots__=['_file']

    def __init__(self, mol2file):
        self._file = mol2file
        self._sasa = None

    def mol2_data(self):
        with open(self._file, 'r') as f:
            yield from f
            # mol2_data = f.readlines()
        # return mol2_data

    def atoms(self):
        atom_lines = []
        collect_atoms = False
        i = 1
        for line in self.mol2_data():
            if '@<TRIPOS>ATOM' in line:
                collect_atoms = True
                continue
            if collect_atoms:
                if line.startswith('@<TRIPOS>'):
                    break
                if line[49] != 'H':
                    if str(line[64:71]).strip()[:3] != 'HOH':
                        # aNum = int(line[0:6])
                        aNum = i
                        aType = str(line[7:14]).strip()
                        aX = float(line[17:26])
                        aY = float(line[28:37])
                        aZ = float(line[40:48])
                        aType2 = str(line[49:56]).strip()
                        rNum = int(line[58:63])
                        rType = str(line[64:71]).strip()
                        aCharge = float(line[72:79])
                        atom_lines.append([aNum, aType, aX, aY, aZ, aType2, rNum, rType, aCharge])
                        i+=1
        return atom_lines
    
    def numAtoms(self):
        return len(self.atoms())
    
    def getNeighbors(self, atom, k=16):
        allNeighbors = []
        for i in self.atoms():
            allNeighbors.append((i, self.calculateDistance(i, atom)))
        sorted_list = sorted(allNeighbors, key=lambda x: x[1])
        return [x for x in sorted_list[:k]]

    def calculateDistance(self, a1, a2):
        return math.sqrt((a1[2]-a2[2])**2+(a1[3]-a2[3])**2+(a1[4]-a2[4])**2)
    
    def adjacencyMatrix(self, atom):
        kNeighbors = self.getNeighbors(atom)
        matrix = np.zeros((len(kNeighbors), len(kNeighbors)))
        for i in range(len(kNeighbors)):
            for j in range(len(kNeighbors)):
                if kNeighbors[i] != kNeighbors[j]:
                    if self.calculateDistance(kNeighbors[i][0], kNeighbors[j][0]) <= 1.55:
                        matrix[i,j] = 1
        return matrix
    
    def featureMatrix(self, atom):
        sasa = self.sasaList()
        data = self.getNeighbors(atom)
        # eDensity = self.electronDensity(data)
        ljP = self.lj_potential(data)
        itsNeighbors = []
        for atomFeature, potential in zip(data, ljP):
            distance = atomFeature[1]
            atomFeature = atomFeature[0]

            atomFeature.append(sasa[atomFeature[0]-1])
            atomFeature[7] = str(atomFeature[7])[0:3]

            for i in potential:
                atomFeature.append(i)
            # del atomFeature[0]
            # del atomFeature[]
            itsNeighbors.append(atomFeature)
        # itsNeighbors = np.array(itsNeighbors)
        return itsNeighbors

    def sasa(self):
        sasa = ShrakeRupley()
        return sasa.compute(np.array(self.atoms()))
    
    def sasaList(self):
        if self._sasa is None:
            sasa = np.array([])
            for i in self.sasa():
                sasa = np.append(sasa, i)
            self._sasa = sasa
        return self._sasa
    

    def predict_secondary_structure(self, atom_coords, atom_type):
        # Calculate distances between atoms
        dist_CA_C = np.linalg.norm(atom_coords['CA'] - atom_coords['C'])
        dist_CA_N = np.linalg.norm(atom_coords['CA'] - atom_coords['N'])
        dist_CA_O = np.linalg.norm(atom_coords['CA'] - atom_coords['O'])
        dist_N_O = np.linalg.norm(atom_coords['N'] - atom_coords['O'])
        
        # Calculate dihedral angles
        dihedral_phi = np.arctan2(np.cross(atom_coords['C'] - atom_coords['CA'], atom_coords['N'] - atom_coords['CA']).dot(np.cross(atom_coords['CA'] - atom_coords['C_prev'], atom_coords['N'] - atom_coords['C'])), np.cross(atom_coords['C'] - atom_coords['CA'], atom_coords['CA'] - atom_coords['N']).dot(np.cross(atom_coords['CA'] - atom_coords['C_prev'], atom_coords['N'] - atom_coords['C'])))
        dihedral_psi = np.arctan2(np.cross(atom_coords['N'] - atom_coords['CA'], atom_coords['C'] - atom_coords['CA']).dot(np.cross(atom_coords['CA'] - atom_coords['N'], atom_coords['C_next'] - atom_coords['CA'])), np.cross(atom_coords['N'] - atom_coords['CA'], atom_coords['CA'] - atom_coords['C']).dot(np.cross(atom_coords['CA'] - atom_coords['N'], atom_coords['C_next'] - atom_coords['CA'])))
        
        # Classify residue as helix, sheet, or loop
        if dist_CA_N < 1.4 and dist_N_O < 1.4 and dihedral_phi < -np.pi/3 and dihedral_psi > np.pi/3:
            return 'H'  # Helix
        elif dist_CA_C < 1.4 and dist_CA_O < 1.4 and dihedral_phi > np.pi/3 and dihedral_psi < -np.pi/3:
            return 'E'  # Sheet
        else:
            return 'C'  # Coil/loop

    # def electronDensity(self, data):
    #     dataStr = ''
    #     for i in data:
    #         goodData = [i[0][5].split('.')[0], str(i[0][2]), str(i[0][3]), str(i[0][4])]
    #         print(goodData)
    #         dataStr += '  '.join(goodData) + '\n'
    #     print(dataStr)
    #     mol = pyscf.gto.M(atom=dataStr, basis='sto-3g')


    #     # Compute the molecular orbitals
    #     mf = pyscf.scf.RHF(mol)
    #     mf.kernel()

    #     # Compute the electron density matrix
    #     dm = mf.make_rdm1()

    #     # Compute the electron density of each atom
    #     electron_density = mol.atom_charges()[:, None] - dm.diagonal()

    #     for i in electron_density:
    #         print(len(i), i)

    #     return data

    import numpy as np

    def lj_potential(self, data):
        """
        Compute the Lennard-Jones potential energy.
        Returns a Matrix with the energies computed between all the neighbors atoms.
        """
        dataMatrix = []
        for i in data:
            dataRow = np.array([i[0][5].split('.')[0], i[0][2], i[0][3], i[0][4]])
            dataMatrix.append(dataRow)
        dataMatrix = np.array(dataMatrix)
        # print(dataMatrix)
        m, n = dataMatrix.shape
        ljP = np.empty([m,m])
        bond = ''
        for i in range(m):
            for j in range(m):
                r = math.sqrt((float(dataMatrix[i, 1])-float(dataMatrix[j, 1]))**2+(float(dataMatrix[i, 2])-float(dataMatrix[j, 2]))**2+(float(dataMatrix[i, 3])-float(dataMatrix[j, 3]))**2)
                if len(str(dataMatrix[i, 0])) == 1:
                    bond = dataMatrix[i, 0]+" -"+dataMatrix[j, 0]
                elif len(str(dataMatrix[i, 0])) == 2:
                    bond = dataMatrix[i, 0]+"-"+dataMatrix[j, 0]
                # print('Bond', bond)
                # with open('parm99.dat') as ff:
                eps = 0
                sig = 0
                # for line in ff:
                for line in self.readAtomES():
                    if line.startswith(bond):
                        eps = float(line[7:12])
                        sig = float(line[16:22])
                        # print(eps, sig)
                        break
                if eps == 0 and sig == 0:
                    if len(str(dataMatrix[i, 0])) == 1:
                        bond = dataMatrix[j, 0]+" -"+dataMatrix[i, 0]
                    elif len(str(dataMatrix[i, 0])) == 2:
                        bond = dataMatrix[j, 0]+"-"+dataMatrix[i, 0]
                    # with open('parm99.dat') as ff:
                        # eps = 0
                        # sig = 0
                        # for line in ff:
                    for line in self.readAtomES():
                        # print('Second:', line)
                        if line.startswith(bond):
                            eps = float(line[7:12])
                            sig = float(line[16:22])
                            # print(eps, sig)
                            break
                if r != 0:
                    ljP[i, j] = "{:.3}".format(4 * eps * ((sig / r) ** 12 - (sig / r) ** 6))
                else:
                    ljP[i, j] = 0
        # print(ljP)
        return ljP
    
    def readAtomES(self):
        with open('parm99.dat') as ff:
            yield from ff



class Mol2ligand():
    def __init__(self, mol2ligand):
        # self.__mol = mol
        self.__ligand = mol2ligand
    def readSolution(self):
        with open(self.__ligand, 'r') as f:
            ligand_data = f.readlines()
        return ligand_data
    def points(self):
        points = []
        collect_atoms = False
        for line in self.readSolution():
            if '@<TRIPOS>ATOM' in line:
                collect_atoms = True
                continue
            if collect_atoms:
                if line.startswith('@<TRIPOS>'):
                    break
                points.append([float(line[18:26]), float(line[29:37]), float(line[40:48])])
            return points
    def SolutionsFeatureMatrix(self, featureMatrix):
        featMatrix = []
        for atom in featureMatrix:
            ligand = 0
            for cavity in self.points():
                dist = ((atom[2]-cavity[0])**2+(atom[3]-cavity[1])**2+(atom[4]-cavity[2])**2)
                if dist < 2:
                    ligand = 1
                    break
            atom.append(ligand)
            featMatrix.append(atom)
        return featMatrix



ATOMIC_RADII = collections.defaultdict(lambda: 2.0)
ATOMIC_RADII.update(
    {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }
)

class ShrakeRupley:
    """Calculates SASAs using the Shrake-Rupley algorithm."""

    def __init__(self, probe_radius=1.40, n_points=960, radii_dict=None):
        """Initialize the class.

        :param probe_radius: radius of the probe in A. Default is 1.40, roughly
            the radius of a water molecule.
        :type probe_radius: float

        :param n_points: resolution of the surface of each atom. Default is 100.
            A higher number of points results in more precise measurements, but
            slows down the calculation.
        :type n_points: int

        :param radii_dict: user-provided dictionary of atomic radii to use in
            the calculation. Values will replace/complement those in the
            default ATOMIC_RADII dictionary.
        :type radii_dict: dict

        >>> sr = ShrakeRupley()
        >>> sr = ShrakeRupley(n_points=960)
        >>> sr = ShrakeRupley(radii_dict={"O": 3.1415})
        """
        if probe_radius <= 0.0:
            raise ValueError(
                f"Probe radius must be a positive number: {probe_radius} <= 0"
            )
        self.probe_radius = float(probe_radius)

        if n_points < 1:
            raise ValueError(
                f"Number of sphere points must be larger than 1: {n_points}"
            )
        self.n_points = n_points

        # Update radii list with user provided lists.
        self.radii_dict = ATOMIC_RADII.copy()
        if radii_dict is not None:
            self.radii_dict.update(radii_dict)

        # Pre-compute reference sphere
        self._sphere = self._compute_sphere()

    def _compute_sphere(self):
        """Return the 3D coordinates of n points on a sphere.

        Uses the golden spiral algorithm to place points 'evenly' on the sphere
        surface. We compute this once and then move the sphere to the centroid
        of each atom as we compute the ASAs.
        """
        n = self.n_points

        dl = np.pi * (3 - 5**0.5)
        dz = 2.0 / n

        longitude = 0
        z = 1 - dz / 2

        coords = np.zeros((n, 3), dtype=np.float32)
        for k in range(n):
            r = (1 - z * z) ** 0.5
            coords[k, 0] = math.cos(longitude) * r
            coords[k, 1] = math.sin(longitude) * r
            coords[k, 2] = z
            z -= dz
            longitude += dl

        return coords

    def compute(self, atoms):
        """Calculate surface accessibility surface area for an entity.

        The resulting atomic surface accessibility values are attached to the
        .sasa attribute of each entity (or atom), depending on the level. For
        example, if level="R", all residues will have a .sasa attribute. Atoms
        will always be assigned a .sasa attribute with their individual values.

        :param entity: input entity.
        """
        # Get atoms and coords
        n_atoms = len(atoms)
        coords = np.array([a[2:5] for a in atoms], dtype=np.float64)

        # Pre-compute atom neighbors using KDTree
        kdt = KDTree(coords, 10)

        # Pre-compute radius * probe table
        radii_dict = self.radii_dict
        radii = np.array([radii_dict[str(a[5].split('.')[0]).upper()] for a in atoms], dtype=np.float64)
        radii += self.probe_radius
        twice_maxradii = np.max(radii) * 2

        # Calculate ASAa
        asa_array = np.zeros((n_atoms, 1), dtype=np.int64)
        ptset = set(range(self.n_points))
        for i in range(n_atoms):
            r_i = radii[i]

            # Move sphere to atom
            s_on_i = (np.array(self._sphere, copy=True) * r_i) + coords[i]
            available_set = ptset.copy()

            # KDtree for sphere points
            kdt_sphere = KDTree(s_on_i, 10)

            # Iterate over neighbors of atom i
            for jj in kdt.search(coords[i], twice_maxradii):
                j = jj.index
                if i == j:
                    continue

                if jj.radius < (r_i + radii[j]):
                    # Remove overlapping points on sphere from available set
                    available_set -= {
                        pt.index for pt in kdt_sphere.search(coords[j], radii[j])
                    }

            asa_array[i] = len(available_set)  # update counts

        # Convert accessible point count to surface area in A**2
        f = radii * radii * (4 * np.pi / self.n_points)
        asa_array = asa_array * f[:, np.newaxis]

        # Set atom .sasa
        sasa = []
        for i, atom in enumerate(atoms):
            sasa.append(asa_array[i, 0])
        return sasa







