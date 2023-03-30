import proteinclass
import numpy as np
import math
import itertools
import sys
from numba import njit
from time import time
from scipy.spatial import cKDTree
from small_sasa import ShrakeRupley2

mol = proteinclass.readMod2('scPDB/1iki_1/protein.mol2')

radii_dict = {
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


def get_coord_in_list(molecule):
    coords = []
    radii = []
    for i in molecule:
        coords.append([float(i[2]), float(i[3]), float(i[4])])
        radii.append(radii_dict[i[5].split('.')[0]])
    return np.array(coords), np.array(radii)

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


def get_angles_and_atoms_3(outside_points, angles_combination, radii, tree_coords):

    bs = open('binding_points.pdb', 'a')
    n = 1

    #radii_inside = radii[tree_coords.query(outside_points)[1]]
    # print(tree_coords.query(outside_points)[0], radii[tree_coords.query(outside_points)[1]]*2)
    mask = tree_coords.query(outside_points)[0] <= radii[tree_coords.query(outside_points)[1]]
    # print(mask)
    for ijk, m in zip(outside_points, mask):
        # print(ijk)
        if m:
            # print('Inside', ijk)
            continue


        theta_collisions = []
        phi_collisions = []
        # points = np.array([np.add(np.array(ijk), np.array(spherical_to_cartesian(d, t, p))) for t, p in angles_combination for d in range(1, 10)])
        for angle in angles_combination:
            t = angle[0]
            p = angle[1]
            points = np.array([np.add(np.array(ijk), np.array(spherical_to_cartesian(d, t, p))) for d in range(2, 15)])
            collisions = tree_coords.query(points)[0] <= radii[tree_coords.query(ijk)[1]]
            if np.any(collisions):
                theta_collisions.append(t)
                phi_collisions.append(p)
        # print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)

        for start_angle in range(0, 181, 30):
            end_angle = start_angle + 181
            range_list = list(range(start_angle, end_angle, 30))
            if all(elem in theta_collisions for elem in range_list):
                for start_angle2 in range(0, 121, 30):
                    end_angle2 = start_angle2 + 241
                    range_list2 = list(range(start_angle2, end_angle2, 30))
                    if all(elem in phi_collisions for elem in range_list2):
                        ijk = np.round(ijk, decimals=3)
                        print(ijk)
                        print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)
                        bs.write("ATOM"+" "*(7-len(str(n)))+str(n)+"  O   HOH X"+" "*(4-len(str(n)))+str(n)+" "*(8-len(str(str(ijk[0]).split('.')[0])))+str(format(ijk[0], ".3f"))+" "*(4-len(str(str(ijk[1]).split('.')[0])))+str(format(ijk[1], ".3f"))+" "*(4-len(str(str(ijk[2]).split('.')[0])))+str(format(ijk[2], ".3f"))+"  1.00 30.00           O  \n")
                        n += 1
                        break
                        #print('yes')
                        
    bs.write("END                                                                             \n")    




coords, radii = get_coord_in_list(mol.atoms())
radii = np.array(radii)

coordsMax = [math.ceil(coords[:,0].max()), math.ceil(coords[:,1].max()), math.ceil(coords[:,2].max())]
coordsMin = [math.ceil(coords[:,0].min()-1), math.ceil(coords[:,1].min()-1), math.ceil(coords[:,2].min()-1)]

print(coordsMax, coordsMin, np.subtract(coordsMax,coordsMin))





xp = np.arange(coordsMin[0], coordsMax[0]+1)
yp = np.arange(coordsMin[1], coordsMax[1]+1)
zp = np.arange(coordsMin[2], coordsMax[2]+1)
grid_points = np.array(list(itertools.product(xp,yp,zp))) # Generate all the grid points






print("starting collision checking")
start = time()


theta = list(range(0, 361, 30))
phi = list(range(0, 361, 30))
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

sasa = ShrakeRupley2()
sasa_values = np.array(sasa.compute(coords, radii))

exposed_coords = []
exposed_radii = []
for s, c, r, sas in zip(sasa_values, coords, radii, sasa_values):
    if sas >= 0:
        exposed_coords.append(c.tolist())
        exposed_radii.append(r)

x = [-2,0,2]
grid_sasa_points = np.array(list(itertools.product(x,x,x))) # Generate all the grid points
surface_grid = []
for ec in exposed_coords:
    for gsp in grid_sasa_points:
        surface_grid.append(np.round(np.add(np.array(ec), gsp), decimals=0))

surface_grid = np.unique(np.array(surface_grid), axis=0)
tree = cKDTree(coords)
get_angles_and_atoms_3(surface_grid, angles_combination, radii, tree)


print(f"finished after {round(time() - start, 2)} seconds")


