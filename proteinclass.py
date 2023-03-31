import math
import numpy as np
import collections
from ctypes import *
from Bio.PDB.kdtrees import KDTree
import pyscf
import chemistry as chem
from scipy.spatial import cKDTree
import itertools
from numba import njit


class readprotein():
   # __slots__=['_file']


    def __init__(self, file):
        self._file = file
        self._sasa = None
        self._mol_dict = None
        self._secondary = None
        self._pocket_tree = None


    def data(self):
        with open(self._file, 'r') as f:
            yield from f
    
    def numAtoms(self):
        return len(self.atoms())
    
    def numAA(self):
        num = 0
        for i in self.atoms():
            if i[6] != num:
                num += 1
        return num
    
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
    

    
    def geometric_binding(self):
        @njit
        def spherical_to_cartesian(r, theta, phi):
            theta = theta/180*math.pi
            phi = phi/180*math.pi
            x = r * math.sin(theta) * math.cos(phi)
            y = r * math.sin(theta) * math.sin(phi)
            z = r * math.cos(theta)
            x = round(x, 3)
            y = round(y, 3)
            z = round(z, 3)
            return np.array([x, y, z])

        coords = []
        radii = []
        for i in self.atoms():
            coords.append([float(i[2]), float(i[3]), float(i[4])])
            # radii.append(radii_dict[i[5].split('.')[0]])
            radii.append(ATOMIC_RADII[i[5].split('.')[0].upper()])
        coords = np.array(coords)
        radii = np.array(radii)
        # return coords, radii
        # return self.data()

        theta = list(range(0, 361, 45))
        phi = list(range(0, 361, 45))
        angles_combination = np.array(list(itertools.product(theta, phi)))
        # print(angles_combination)
        visited_angles = []
        unique_angles = []
        for angle in angles_combination:
            t = angle[0]
            p = angle[1]
            #### Do not repeat same direction vectors ####
            if spherical_to_cartesian(1, t, p).tolist() not in visited_angles:
                visited_angles.append(spherical_to_cartesian(1, t, p).tolist())
                unique_angles.append([t, p])
            else:
                pass
            #### END BLOCK ####
        angles_combination = unique_angles

        sasa_values = self.sasaList()
        # print(sasa)

        # sasa = ShrakeRupley2()
        # sasa_values = np.array(sasa.compute(coords, radii))
        exposed_coords = []
        exposed_radii = []
        for c, r, sas in zip(coords, radii, sasa_values):
            if sas >= 0:
                exposed_coords.append(c.tolist())
                exposed_radii.append(r)

        x = [-3,0,3]
        grid_sasa_points = np.array(list(itertools.product(x,x,x))) # Generate all the grid points
        surface_grid = []
        for ec in exposed_coords:
            for gsp in grid_sasa_points:
                surface_grid.append(np.round(np.add(np.array(ec), gsp), decimals=0))

        surface_grid = np.unique(np.array(surface_grid), axis=0)
        tree = cKDTree(coords)

        pocket_points = []
        n = 1
        mask = tree.query(surface_grid)[0] <= radii[tree.query(surface_grid)[1]]

        print(mask)

        surface_grid = surface_grid[np.where(~mask)]

        print(surface_grid)

        # for ijk, m in zip(surface_grid, mask):
        #     if m:
        #         continue
        #     theta_collisions = []
        #     phi_collisions = []
        #     for angle in angles_combination:
        #         t = angle[0]
        #         p = angle[1]
        #         points = np.array([np.add(np.array(ijk), np.array(spherical_to_cartesian(d, t, p))) for d in range(2, 15)])
        #         collisions = tree.query(points)[0] <= radii[tree.query(ijk)[1]]
        #         if np.any(collisions):
        #             theta_collisions.append(t)
        #             phi_collisions.append(p)
        #     for start_angle in range(0, 181, 30):
        #         end_angle = start_angle + 181
        #         range_list = list(range(start_angle, end_angle, 30))
        #         if all(elem in theta_collisions for elem in range_list):
        #             for start_angle2 in range(0, 121, 30):
        #                 end_angle2 = start_angle2 + 241
        #                 range_list2 = list(range(start_angle2, end_angle2, 30))
        #                 if all(elem in phi_collisions for elem in range_list2):
        #                     ijk = np.round(ijk, decimals=3)
        #                     pocket_points.append(ijk)
        #                     print(ijk, start_angle, start_angle2)
        #                     n += 1
        #                     break
        


        # Compute queries for all points
        # dists, indices = tree.query(np.array([np.add(np.array(surface_grid), np.array(spherical_to_cartesian(d, t, p))) for d in range(4, 15) for t, p in angles_combination]))

        # for i, (ijk, m) in enumerate(zip(surface_grid, mask)):
        #     if m:
        #         continue
        #     theta_collisions = set()
        #     phi_collisions = set()
        #     for j, (t, p) in enumerate(angles_combination):
        #         collisions = dists[(j * 13):((j + 1) * 13)] <= radii[indices[(j * 13):((j + 1) * 13)]]
        #         if np.any(collisions):
        #             theta_collisions.add(t)
        #             phi_collisions.add(p)
        #     for start_angle in range(0, 181, 30):
        #         end_angle = start_angle + 181
        #         range_set = set(range(start_angle, end_angle, 30))
        #         if range_set.issubset(theta_collisions):
        #             for start_angle2 in range(0, 121, 30):
        #                 end_angle2 = start_angle2 + 241
        #                 range_set2 = set(range(start_angle2, end_angle2, 30))
        #                 if range_set2.issubset(phi_collisions):
        #                     ijk = np.round(ijk, decimals=3)
        #                     pocket_points.append(ijk)
        #                     print(ijk, start_angle, start_angle2)
        #                     n += 1
        #                     break






        return cKDTree(pocket_points)
    
    def pocketTree(self):
        if self._pocket_tree is None:
            self._pocket_tree = self.geometric_binding()
        return self._pocket_tree
    
    def featureMatrix(self, atom):
        '''
        Get the feature Matrix for the atom and its neighbors.
        ---
        Position and description:
        []


        '''

        pocket_tree = self.pocketTree()

        print(pocket_tree)

        sasa = self.sasaList()
        for i in range(0, len(sasa)):
            sasa[i] = (sasa[i]-np.min(sasa))/(np.max(sasa)-np.min(sasa))
        data = self.getNeighbors(atom)
        # print(data)
        directions = self.directions(data)
        # eDensity = self.electronDensity(data)
        ljP = self.lj_potential(data)
        secondary = self.secondaryList()
        sec_number = []
        # print('-'*10+'features'+'-'*10)
        for i in secondary:
            if i == '?': # Loop
                sec_number.append('0')
            if i == 'E': # Beta sheet
                sec_number.append('1')
            if i == 'H': # Alpha helix
                sec_number.append('2')
        featMat = []
        for atomFeature, potential, direction in zip(data, ljP, directions):
            distance = atomFeature[1]
            atomFeature = atomFeature[0]


            atomFeature.append(sasa[atomFeature[0]-1])


            for i in potential:
                atomFeature.append(i)


            for d in [chem.dictionary_kd_hydrophobicity, chem.dictionary_ww_hydrophobicity, chem.dictionary_hh_hydrophobicity, chem.dictionary_mf_hydrophobicity, chem.dictionary_tt_hydrophobicity]:
                if atomFeature[7] not in d:
                    for i in range(5):
                        atomFeature.append(0)
                else:
                    atomFeature.append(chem.dictionary_kd_hydrophobicity[atomFeature[7]])
                    atomFeature.append(chem.dictionary_ww_hydrophobicity[atomFeature[7]])
                    atomFeature.append(chem.dictionary_hh_hydrophobicity[atomFeature[7]])
                    atomFeature.append(chem.dictionary_mf_hydrophobicity[atomFeature[7]])
                    atomFeature.append(chem.dictionary_tt_hydrophobicity[atomFeature[7]])
            
            for i in direction:
                atomFeature.append(i)
            

            # print(sec_number)
            # print(atomFeature[6]-1, len(sec_number))
            # print(sec_number[atomFeature[6]-1])
            atomFeature.append(sec_number[atomFeature[6]-1])


            ### DELETE ALL UNDESIRED ELEMENTS ###
            del atomFeature[0:8]



            ### CONVERT ALL ELEMENTS TO FLOAT ###
            for i in  range(0, len(atomFeature)):
                atomFeature[i] = float(atomFeature[i])



            # del atomFeature[0]
            # del atomFeature[0]
            # del atomFeature[3]
            # del atomFeature[3]
            # del atomFeature[3]
            # # print(atomFeature)
            # # del atomFeature[]


            # del atomFeature[0:3]
            # del atomFeature[2:]




            # print(atomFeature)


            featMat.append(atomFeature)

        print(featMat)
        return featMat


    def directions(self, data):
        given_atom = data[0][0][2:5]
        neighbor_atoms = []
        for i in data:
            neighbor_atoms.append(i[0][2:5])
        given_atom = np.array(given_atom)
        neighbor_atoms = np.array(neighbor_atoms)
        directions = neighbor_atoms-given_atom


        return directions


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
    
    def coords(self):
        if self._mol_dict is None:
            moleculeList = []
            c = []
            smallDict = {}
            for i in self.atoms():
                if i[6] not in c:
                    if i[1] == 'N':
                        smallDict['N'] = np.array([i[2], i[3], i[4]])
                    elif i[1] == 'CA':
                        smallDict['CA'] = np.array([i[2], i[3], i[4]])
                    elif i[1] == 'C':
                        smallDict['C'] = np.array([i[2], i[3], i[4]])
                    elif i[1] == 'O':
                        smallDict['O'] = np.array([i[2], i[3], i[4]])
                if len(smallDict) == 4:
                    smallDict['R'] =i[6]
                    moleculeList.append(smallDict)
                    smallDict = {}
                    c.append(i[6])
            # print(c, len(c))
            self._mol_dict = moleculeList
        return self._mol_dict


    def secondaryList(self):
        if self._secondary is None:
            sec = np.array([])
            for i in self.calculate_secondary_structure(self.coords()):
                sec = np.append(sec, i)
            while self.numAA() > len(sec):
                sec = np.append(sec, '?')
            # print(self.numAA())
            # print(len(sec))
            self._secondary = sec
        return self._secondary


    def calculate_secondary_structure(self, coords):
        # Define the hydrogen bonding distance cutoffs for each type of interaction
        hbond_cutoffs = {'helix': 3.4, 'sheet': 3.2}
        
        # Define the angles for each type of secondary structure
        helix_angles = {'phi': 92, 'psi': 98}
        sheet_angles = {'C-alpha-C': 139, 'C-N-C-alpha': 135}
        
        # Initialize the secondary structure assignment for each residue
        sec_struct = ['?' for i in range(len(coords))]
        
        # Iterate over each residue and compare its hydrogen bonding distance to adjacent residues
        for i in range(len(coords)-3):
            for j in range(i+1, len(coords)-3):
                if (np.linalg.norm(coords[i]['O'] - coords[j]['N']) < hbond_cutoffs['helix'] and (abs(coords[j]['R']-coords[i]['R']) >= 3 and abs(coords[j]['R']-coords[i]['R']) <= 5)) or (np.linalg.norm(coords[i]['N'] - coords[j]['O']) < hbond_cutoffs['helix'] and (abs(coords[j]['R']-coords[i]['R']) >= 3 and abs(coords[j]['R']-coords[i]['R']) <= 5)):
                    # print(coords[i]['R']+2, coords[j]['R']+2)
                    phi1 = self.calculate_angle(coords[i]['C'], coords[i]['N'], coords[i+1]['CA'], coords[i+1]['C'])
                    psi1 = self.calculate_angle(coords[i]['N'], coords[i]['CA'], coords[i+1]['C'], coords[i+1]['N'])
                    # print(phi1, psi1)
                    phi2 = self.calculate_angle(coords[i+1]['C'], coords[i+1]['N'], coords[i+2]['CA'], coords[i+2]['C'])
                    psi2 = self.calculate_angle(coords[i+1]['N'], coords[i+1]['CA'], coords[i+2]['C'], coords[i+2]['N'])
                    # print(phi2, psi2)
                    phi3 = self.calculate_angle(coords[i+2]['C'], coords[i+2]['N'], coords[i+3]['CA'], coords[i+3]['C'])
                    psi3 = self.calculate_angle(coords[i+2]['N'], coords[i+2]['CA'], coords[i+3]['C'], coords[i+3]['N'])
                    # print(abs((angle1 - helix_angles['C-alpha-C'])/helix_angles['C-alpha-C']), abs((angle2 - helix_angles['C-N-C-alpha'])/helix_angles['C-N-C-alpha']))
                    if abs((phi1 - helix_angles['phi'])/helix_angles['phi']) < 0.35 and abs((psi1 - helix_angles['psi'])/helix_angles['psi']) < 0.35 and abs((phi2 - helix_angles['phi'])/helix_angles['phi']) < 0.35 and abs((psi2 - helix_angles['psi'])/helix_angles['psi']) < 0.35 and abs((phi3 - helix_angles['phi'])/helix_angles['phi']) < 0.35 and abs((psi3 - helix_angles['psi'])/helix_angles['psi']) < 0.35:
                        sec_struct[i] = 'H'
                        sec_struct[i+1] = 'H'
                        sec_struct[i+2] = 'H'
        for i in range(len(coords)-1):
            for j in range(i+1, len(coords)-1):
                if (np.linalg.norm(coords[i]['O'] - coords[j]['N']) < hbond_cutoffs['sheet'] or np.linalg.norm(coords[i]['N'] - coords[j]['O']) < hbond_cutoffs['sheet']) and sec_struct[i+2] != 'H' and sec_struct[i+1] != 'H' and sec_struct[i] != 'H' and sec_struct[j] != 'H' and coords[j]['R'] != coords[i]['R']+1:
                    # print('E?', coords[i]['R'], coords[j]['R'], np.linalg.norm(coords[i]['O'] - coords[j]['N']), np.linalg.norm(coords[i]['N'] - coords[j]['O']), sec_struct[i], sec_struct[j])
                    if np.linalg.norm(coords[i]['O'] - coords[j]['N']) < np.linalg.norm(coords[i]['N'] - coords[j]['O']):
                        cac = self.midpoint(coords[j-1]['C'], coords[j]['CA']) # Midpoint between CA and C
                        hydrogen = self.placeH(coords[j]['N'], cac) # Place an hydrogen (prefectly geometric)
                        angle = self.Hangle(coords[i]['O'], hydrogen, coords[j]['N'])
                    else:
                        cac = self.midpoint(coords[i-1]['C'], coords[i]['CA']) # Midpoint between CA and C
                        hydrogen = self.placeH(coords[i]['N'], cac) # Place an hydrogen (prefectly geometric)
                        angle = self.Hangle(coords[i]['N'], hydrogen, coords[j]['O'])
                    # print('angle', angle)
                    # Check if the angle between the C-alpha, C, and N atoms is within the sheet range
                    phi = self.calculate_angle(coords[i]['C'], coords[i]['N'], coords[i+1]['CA'], coords[i+1]['C'])
                    psi = self.calculate_angle(coords[i]['N'], coords[i]['CA'], coords[i+1]['C'], coords[i+1]['N'])
                    # print(phi, psi)
                    # if abs((phi - sheet_angles['C-alpha-C'])/sheet_angles['C-alpha-C']) < 0.35 and abs((psi - sheet_angles['C-N-C-alpha'])/sheet_angles['C-N-C-alpha']) < 0.35:
                    if abs((angle - 180)/180) < 0.15:
                        sec_struct[i] = 'E'
                        sec_struct[j] = 'E'
        for i in range(1, len(sec_struct)-1):
            if sec_struct[i] == '?' and sec_struct[i-1] == 'E' and sec_struct[i+1] == 'E':
                sec_struct[i] = 'E'
        # c = 1
        # for i in sec_struct:
        #     if i == 'E':
        #         print('E', c)
        #     if i == 'H':
        #         print('H', c)
        #     c+=1
        
        return sec_struct
    
    def placeH(self, n, mp):
        direction = n-mp
        h = n+direction
        return h
    
    def midpoint(self, p1, p2):
        return (p1+p2)/2

    def Hangle(self, p0, p1, p2):
        ba = p0 - p1
        bc = p2 - p1
        cos = np.dot(ba, bc)/(np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cos)
        angle = np.degrees(angle)

        return angle

    def calculate_angle(self, p0, p1, p2, p3):
        b0 = p0 - p1
        b1 = p2 - p1
        b2 = p3 - p2

        # https://leimao.github.io/blog/Dihedral-Angles/

        cos = np.dot(np.cross(b0, b1), np.cross(b1, b2))/(np.linalg.norm(np.cross(b0, b1)) * np.linalg.norm(np.cross(b1, b2)))
        sin = np.dot(np.cross(np.cross(b0, b1), np.cross(b1, b2)), b1)/(np.linalg.norm(np.cross(b0, b1)) * np.linalg.norm(np.cross(b1, b2)) * np.linalg.norm(b1))

        angle = np.arctan2(sin, cos)
        angle = np.rad2deg(angle)

        return angle


    # def calculate_angle(self, p1, p2, p3):
    #     v1 = p1 - p2
    #     v2 = p3 - p2
    #     cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    #     theta = np.arccos(cos_theta)
    #     return theta * 180 / np.pi






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




class readMod2(readprotein):
    '''
    Explicit to read the mol2 file
    '''


    def atoms(self):
        atom_lines = []
        collect_atoms = False
        i = 1
        for line in self.data():
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
                        atom_lines.append([aNum, aType, aX, aY, aZ, aType2, rNum, rType[0:3], aCharge])
                        i+=1
        return atom_lines




class readpdb(readprotein):
    '''
    Explicit to read the pdb file
    '''
    pass


class Mol2ligand():


    def __init__(self, mol2ligand):
        self._ligand = mol2ligand


    def readSolution(self):
        with open(self._ligand, 'r') as fl:
            yield from fl
            # ligand_data = fl.readlines()
        # return ligand_data
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
                points.append([float(line[17:25]), float(line[27:35]), float(line[37:45])])
        return points
        
    def SolutionsFeatureMatrix(self, matrix, sasa):
        solList = []
        for atom in matrix:
            atom = atom[0]
            # print(atom)
            ligand = 0
            for cavity in self.points():
                dist = math.sqrt((atom[2]-cavity[0])**2+(atom[3]-cavity[1])**2+(atom[4]-cavity[2])**2)
                # print([atom[2],cavity[0],atom[3],cavity[1],atom[4],cavity[2]])
                # print(dist)

                if dist < 3.5 and sasa[atom[0]-1] > 0:
                    ligand = 1
                    break
            solList.append(ligand)
        return solList






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

class ShrakeRupley2:
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


    def compute(self, coords, radii):
        """Calculate surface accessibility surface area for an entity.


        The resulting atomic surface accessibility values are attached to the
        .sasa attribute of each entity (or atom), depending on the level. For
        example, if level="R", all residues will have a .sasa attribute. Atoms
        will always be assigned a .sasa attribute with their individual values.


        :param entity: input entity.
        """
        # Get atoms and coords
        n_atoms = len(coords)
        # coords = np.array([a[2:5] for a in atoms], dtype=np.float64)


        # Pre-compute atom neighbors using KDTree
        kdt = KDTree(coords, 10)


        # Pre-compute radius * probe table
        # radii_dict = self.radii_dict
        # radii = np.array([radii_dict[str(a[5].split('.')[0]).upper()] for a in atoms], dtype=np.float64)
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
        for i, atom in enumerate(coords):
            sasa.append(asa_array[i, 0])
        return sasa