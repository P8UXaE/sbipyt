import proteinclass
import numpy as np
import math
import itertools
import sys
from numba import njit
from time import time
from scipy.spatial import cKDTree
# from angles import get_angles_and_atoms_3
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

# @njit
# def check_inside(grid_points, coords, radii):
#     outsideList = []
#     for g in grid_points:
#         outside = True
#         # print(g)
#         coords_reduced = coords[np.abs(coords[:, 0] - g[0]) <= 2] # Take the coordinates where the distance x-x is 2 or less, because if it is grater, it won't be inside any atom.
#         for c, r in zip(coords_reduced, radii):
#             if np.linalg.norm(g-c) <= r:
#                 # print(g, np.linalg.norm(g-c))
#                 outside = False
#                 break
#         if outside:
#             outsideList.append(g)
#     return outsideList



# def check_inside_kdtree(grid_points, coords, radii):
#     tree = cKDTree(coords)
#     outsideList = []
#     for g in grid_points:
#         outside = True
#         # print(g)
#         tree.query(g)
#         if tree.query(g)[0] <= radii[tree.query(g)[1]]:
#             # print(g, np.linalg.norm(g-c))
#             outside = False
#             break
#         if outside:
#             outsideList.append(g)
#     return outsideList


# def get_angles_and_atoms_3(outside_points, angles_combination, radii, tree_coords):

#     for combination in outside_points:
#         # i = combination[0]
#         # j = combination[1]
#         # k = combination[2]
#         # print(i, j, k)
#         print(combination)
#         if tree_coords.query(combination)[0] <= radii[tree_coords.query(combination)[1]]: # check if we are inside an atom
#             print('Inside')
#             continue
#         theta_collisions = []
#         phi_collisions = []
#         for angle in angles_combination:
#             t = angle[0]
#             p = angle[1]
#             #collistion_detected = False
#             for d in range(0, 25, 1):
#                 secondary_point = np.add(np.array(combination), np.array(spherical_to_cartesian(d, t, p)))
#                 if tree_coords.query(secondary_point)[0] <= radii[tree_coords.query(combination)[1]]:
#                     theta_collisions.append(t)
#                     phi_collisions.append(p)
#                     break

#         print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)


def get_angles_and_atoms_3(outside_points, angles_combination, radii, tree_coords):
    #radii_inside = radii[tree_coords.query(outside_points)[1]]
    # print(tree_coords.query(outside_points)[0], radii[tree_coords.query(outside_points)[1]]*2)
    mask = tree_coords.query(outside_points)[0] <= radii[tree_coords.query(outside_points)[1]]*1.5
    # print(mask)
    for ijk, m in zip(outside_points, mask):
        # print(ijk)
        if m:
            # print('Inside', ijk)
            continue


        theta_collisions = []
        phi_collisions = []
        points = np.array([np.add(np.array(ijk), np.array(spherical_to_cartesian(d, t, p))) for t, p in angles_combination for d in range(1, 25)])
        for angle in angles_combination:
            t = angle[0]
            p = angle[1]
            points = np.array([np.add(np.array(ijk), np.array(spherical_to_cartesian(d, t, p))) for d in range(1, 25)])
            collisions = tree_coords.query(points)[0] <= radii[tree_coords.query(ijk)[1]]
            if np.any(collisions):
                theta_collisions.append(t)
                phi_collisions.append(p)
        # print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)




        for start_angle in range(0, 181, 45):
            end_angle = start_angle + 181
            range_list = list(range(start_angle, end_angle, 45))
            if all(elem in theta_collisions for elem in range_list):
                for start_angle2 in range(0, 181, 45):
                    end_angle2 = start_angle2 + 181
                    range_list2 = list(range(start_angle2, end_angle2, 45))
                    if all(elem in phi_collisions for elem in range_list2):
                        print(ijk)
                        break
                        #print('yes')
                        # print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)



        # angles = []
        # # print(np.array([np.add(np.array(ijk), np.array([spherical_to_cartesian(d, t, p) for d in range(1, 20, 2)])) for t, p in angles_combination]))
        # mask2 = tree_coords.query(np.array([np.add(np.array(ijk), np.array([spherical_to_cartesian(d, t, p) for d in range(1, 20, 2)])) for t, p in angles_combination])) <= radii[tree_coords.query(ijk)[1]]
        # # print(mask2.flatten().tolist())
        # #sys.exit()
        # # if mask2:
        # angles.append(np.where(mask2.flatten()))
        # print(angles)


        




coords, radii = get_coord_in_list(mol.atoms())
radii = np.array(radii)

coordsMax = [math.ceil(coords[:,0].max()), math.ceil(coords[:,1].max()), math.ceil(coords[:,2].max())]
coordsMin = [math.ceil(coords[:,0].min()-1), math.ceil(coords[:,1].min()-1), math.ceil(coords[:,2].min()-1)]

print(coordsMax, coordsMin, np.subtract(coordsMax,coordsMin))





xp = np.arange(coordsMin[0], coordsMax[0]+1)
yp = np.arange(coordsMin[1], coordsMax[1]+1)
zp = np.arange(coordsMin[2], coordsMax[2]+1)
grid_points = np.array(list(itertools.product(xp,yp,zp))) # Generate all the grid points


# print("starting inside atom point checker")
# start = time()
# # outside_points = check_inside(grid_points, coords, radii) # Get which grid points are not inside an atom.
# outside_points = check_inside_kdtree(grid_points, coords, radii)
# print(f"finished after {round(time() - start, 2)} seconds")



# for g in grid_points:
#     for c, r in zip(coords, radii):
#         if abs(g[0]-c[0])<2 or abs(g[1]-c[1])<2 or abs(g[2]-c[2])<2:
#             if np.linalg.norm(g-c) <= r:
#                 print(g, np.linalg.norm(g-c))
#                 break


# WITH SASA #










print("starting collision checking")
start = time()


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


maxDistanceConsideredBinding = 30

# tree = cKDTree(coords)
# get_angles_and_atoms_3(grid_points, angles_combination, radii, tree)

# print(radii)
# print(coords)

sasa = ShrakeRupley2()
sasa_values = np.array(sasa.compute(coords, radii))

exposed_coords = []
exposed_radii = []
for s, c, r, sas in zip(sasa_values, coords, radii, sasa_values):
    if sas >= 15:
        exposed_coords.append(c.tolist())
        exposed_radii.append(r)

xs = [-3, 0, 3]
ys = [-3, 0, 3]
zs = [-3, 0, 3]
x = [-4,0,4]
grid_sasa_points = np.array(list(itertools.product(x,x,x))) # Generate all the grid points
# print(grid_sasa_points)
surface_grid = []
for ec in exposed_coords:
    for gsp in grid_sasa_points:
        surface_grid.append(np.add(np.array(ec), gsp))

surface_grid = np.array(surface_grid)
tree = cKDTree(coords)
get_angles_and_atoms_3(surface_grid, angles_combination, radii, tree)




# outside_points = [[-13, -32, 14]]

# for combination in outside_points:
#     i = combination[0]
#     j = combination[1]
#     k = combination[2]
#     print(i, j, k)
        
#     theta_collisions = []
#     phi_collisions = []
    
#     visited_angles = []
#     for t in theta:
#         for p in phi: # This iterates over the anlges and directions of the point
#             collistion_detected = False
#             #### Do not repeat same direction vectors ####
#             if np.round(spherical_to_cartesian(1, t, p), decimals = 3).tolist() not in visited_angles:
#                 visited_angles.append(np.round(spherical_to_cartesian(1, t, p), decimals = 3).tolist())
#             else:
#                 break
#             #### END BLOCK ####

#             coords_reduced = coords[np.abs(coords[:, 0] - i) <= 30]
#             radii_reduced = radii[np.abs(coords[:, 0] - i) <= 30]
#             # [coord for coord in coords if abs(coord[0] - i) <= 30]
#             # [radii[i] for i, coord in enumerate(coords) if abs(coord[0] - i) <= 30]
#             for c2, r2 in zip(coords_reduced, radii_reduced): # Iterate over all the atom points and its radii
#                 #### GENERATE FIRST POSITION ####
#                 d = 1 # This is the strarting distance of the sphere.
#                 secondary_point = np.add(np.array([i, j, k]), spherical_to_cartesian(d, t, p))
#                 #print(t, p, np.round(spherical_to_cartesian(d, t, p), decimals = 3).tolist())
#                 #### END BLOCK ####

#                 while coordsMin[0] <= secondary_point[0] <= coordsMax[0] and coordsMin[1] <= secondary_point[1] <= coordsMax[1] and coordsMin[2] <= secondary_point[2] <= coordsMax[2] and d <= maxDistanceConsideredBinding:
#                     if np.linalg.norm(secondary_point - c2) <= r2: # Check if there is any collision between an atom and the points we are throwing
#                         theta_collisions.append(t)
#                         phi_collisions.append(p)
#                         collistion_detected = True
#                         break
#                     d+=1
#                     secondary_point = np.add(np.array([i, j, k]), spherical_to_cartesian(d, t, p))
                    
#                 if collistion_detected:
#                     break
#     print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)
        
# sys.exit()


    
#     print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)
#     for start_angle in range(0, 181, 45):
#         end_angle = start_angle + 180
#         range_list = list(range(start_angle, end_angle, 45))
#         if all(elem in theta_collisions for elem in range_list):
#             print('yes')
#             print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)
#             break



print(f"finished after {round(time() - start, 2)} seconds")


