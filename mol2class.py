import math
import numpy as np
import collections
from ctypes import *
from kdtrees import KDTree
# import sys
# sys.path.append('kdtrees.so')
# from kdtrees import KDTree

class readMol2():

    def __init__(self, mol2file):
        self.__file = mol2file

    def mol2_data(self):
        with open(self.__file, 'r') as f:
            mol2_data = f.readlines()
        return mol2_data

    def atoms(self):
        atom_lines = []
        collect_atoms = False
        for line in self.mol2_data():
            if '@<TRIPOS>ATOM' in line:
                collect_atoms = True
                continue
            if collect_atoms:
                if line.startswith('@<TRIPOS>'):
                    break
                if line[49] != 'H':
                    if str(line[64:71]).strip()[:3] != 'HOH':
                        aNum = int(line[0:6])
                        aType = str(line[7:14]).strip()
                        aX = float(line[17:26])
                        aY = float(line[28:37])
                        aZ = float(line[40:48])
                        aType2 = str(line[49:56]).strip()
                        rNum = int(line[58:63])
                        rType = str(line[64:71]).strip()
                        aCharge = float(line[72:79])
                        atom_lines.append([aNum, aType, aX, aY, aZ, aType2, rNum, rType, aCharge])
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
        sasa = np.array([])
        for i in self.sasa():
            sasa = np.append(sasa, i)
        data = self.getNeighbors(atom)
        itsNeighbors = []
        for atomFeature in data:
            distance = atomFeature[1]
            atomFeature = atomFeature[0]
            

            atomFeature.append(sasa[atomFeature[0]-1])
            atomFeature[7] = str(atomFeature[7])[0:3]


            # del atomFeature[0]
            # del atomFeature[]

            itsNeighbors.append(atomFeature)
        # itsNeighbors = np.array(itsNeighbors)
        return itsNeighbors

    def sasa(self):
        sasa = ShrakeRupley()
        return sasa.compute(np.array(self.atoms()))
    

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







