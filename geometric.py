import proteinclass
import numpy as np
import math
import itertools
import sys
from numba import njit
from time import time
from scipy.spatial import cKDTree
from small_sasa import ShrakeRupley2
from scipy.spatial import Voronoi, Delaunay
from scipy.integrate import nquad

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

def area_of_triangle(vertices):
    """
    Compute the area of a triangle given its three vertices (a, b, c).
    """
    a, b, c = vertices
    ab = b - a
    ac = c - a
    return np.linalg.norm(np.cross(ab, ac)) / 2


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
            collisions = tree_coords.query(points)[0] <= radii[tree_coords.query(points)[1]]
            # print(tree_coords.query(points)[0], radii[tree_coords.query(points)[1]])
            # print(collisions)
            if np.any(collisions):
                theta_collisions.append(t)
                phi_collisions.append(p)
        print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)



        theta = np.deg2rad(theta_collisions)
        phi = np.deg2rad(phi_collisions)

        def cap_area(theta, phi=2*np.pi):
            '''
            theta goes from 0 - π , 0 - 180
            phi goes from 0 - 2π , 0 - 360
            cap_area(np.pi/4, np.pi)
            '''
            a = (2*np.pi*(1-np.cos(theta)))*phi/(2*np.pi)
            return a
        def segment_area(theta1, theta2, phi):
            '''
            phi
            theta2 - theta1, theta1, theta2 != 0, 180
            segment_area(np.pi/4, np.pi/2, np.pi)
            '''
            a2 = (2*np.pi*(1-np.cos(theta2)))
            a1 = (2*np.pi*(1-np.cos(theta1)))
            a = (a2-a1)*phi/(2*np.pi)
            return a

        pairs = np.array([(t, p) for t, p in zip(theta_collisions, phi_collisions)])

        surface = 0
        for th in pairs:
            if 45 <= th[0] <= 90:
                square_opp_corner = np.add(th, np.array([45, 45]))
                if square_opp_corner[1] == 360:
                    square_opp_corner[1] = 0
                if np.any(np.all(pairs == square_opp_corner, axis=1)):
                    surface += segment_area(th[0]*np.pi/180, square_opp_corner[0]*np.pi/180, np.pi/4)
            if abs(th[0]) == 0:
                side_corner = [np.any(np.all(pairs == np.add(th, np.array([45, x])), axis=1)) for x in range(0, 360, 45)]
                for i in range(len(side_corner)-1):
                    if side_corner[i] == True and side_corner[i+1] == True:
                        surface += cap_area(np.pi/4, np.pi/4)
                if side_corner[0] == True and side_corner[7] == True:
                    surface += cap_area(np.pi/4, np.pi/4)
            if abs(th[0]) == 180:
                side_corner = [np.any(np.all(pairs == np.add(th, np.array([-45, x])), axis=1)) for x in range(0, 360, 45)]
                for i in range(len(side_corner)-1):
                    if side_corner[i] == True and side_corner[i+1] == True:
                        surface += cap_area(np.pi/4, np.pi/4)
                if side_corner[0] == True and side_corner[7] == True:
                    surface += cap_area(np.pi/4, np.pi/4)

        # print(surface)

        if surface >= 2*np.pi:
            print("The points cover at least half of the sphere.", surface)
            bs.write("ATOM"+" "*(7-len(str(n)))+str(n)+"  O   HOH X"+" "*(4-len(str(n)))+str(n)+" "*(8-len(str(str(ijk[0]).split('.')[0])))+str(format(ijk[0], ".3f"))+" "*(4-len(str(str(ijk[1]).split('.')[0])))+str(format(ijk[1], ".3f"))+" "*(4-len(str(str(ijk[2]).split('.')[0])))+str(format(ijk[2], ".3f"))+"  1.00 30.00           O  \n")
            n += 1
        else:
            print("The points do not cover at least half of the sphere.")

    bs.write("END                                                                             \n")



        # # Define the set of points on the sphere
        # points = np.array([[t, p] for t, p in zip(theta, phi)])

        # # Convert spherical coordinates to Cartesian coordinates
        # cartesian_points = np.array([spherical_to_cartesian(1, theta, phi) for theta, phi in points])

        # # Compute the Delaunay triangulation of the Cartesian points
        # tri = Delaunay(cartesian_points)

        # # Compute the areas of all the triangles that are covered by at least one point
        # triangle_areas = []
        # for simplex in tri.simplices:
        #     print(simplex)
        #     vertices = cartesian_points[simplex]
        #     print(vertices)
        #     vertex_covered = False
        #     for vertex in vertices:
        #         if np.any((cartesian_points == vertex).all(axis=1)):
        #             vertex_covered = True
        #             break
        #     if vertex_covered:
        #         triangle_area = area_of_triangle(vertices[:3])
        #         triangle_areas.append(triangle_area)

        # # Sum up the areas of all the covered triangles to obtain the total area covered by the points
        # total_area = np.sum(triangle_areas)

        # # Compute the total area of the sphere
        # total_sphere_area = 4 * np.pi

        # # Compute the percentage of the sphere's area covered by the points
        # coverage = total_area / total_sphere_area * 100

        # print("The points cover {:.2f}% of the sphere's area.".format(coverage))


        # def integrand(theta, phi):
        #     return np.sin(theta)

        # # Define the limits of integration
        # theta_limits = [0, np.pi] # integrate over the range [0, pi] for theta
        # phi_limits = [0, 2*np.pi] # integrate over the range [0, 2*pi] for phi

        # # Perform the integration
        # result, _ = nquad(integrand, [phi_limits, theta_limits])

        # # Compute the percentage of the sphere's area
        # sphere_area = 4*np.pi # total area of the sphere
        # portion_area = result / sphere_area * 100

        # print("The portion of the sphere's area is: {:.2f}%".format(portion_area))


        # sum_area = 0

        # # Compute the total surface area of the sphere
        # sphere_area = 4 * np.pi * r**2

        # # Loop over the collision angles and compute the portion of the sphere's area covered by each collision angle
        # for i in range(len(theta_collisions)):
        #     # Extract the polar and azimuthal angles for this collision
        #     theta = np.deg2rad(theta_collisions[i])
        #     phi = np.deg2rad(phi_collisions[i])

        #     # Compute the solid angle for this collision
        #     solid_angle = 2 * np.pi * (1 - np.cos(theta))

        #     # Compute the portion of the sphere's area covered by this collision
        #     portion_area = solid_angle / sphere_area * 100

        #     sum_area += portion_area

        #     # Print the result for this collision
        #     # print("Collision {}: Theta = {} deg, Phi = {} deg, covers {:.2f}% of the sphere's area.".format(i+1, theta_collisions[i], phi_collisions[i], portion_area))

        # print(sum_area)


        # print('theta:', theta, 'Phi:', phi)

        # # Convert theta and phi to Cartesian coordinates
        # x = np.sin(theta) * np.cos(phi)
        # y = np.sin(theta) * np.sin(phi)
        # z = np.cos(theta)





        # # Define the radius of the sphere
        # r = 1.0

        # # Compute the total surface area of the sphere
        # sphere_area = 4 * np.pi * r**2

        # # Initialize the total portion area to zero
        # total_portion_area = 0

        # # Loop over the collision angles and compute the portion of the sphere's area covered by each collision
        # for i in range(len(theta_collisions)):
        #     # Extract the polar and azimuthal angles for this collision
        #     theta = np.deg2rad(theta_collisions[i])
        #     phi = np.deg2rad(phi_collisions[i])

        #     # Compute the solid angle for this collision
        #     solid_angle = 2 * np.pi * (1 - np.cos(theta))

        #     # Compute the portion of the sphere's area covered by this collision
        #     portion_area = solid_angle / sphere_area * 100

        #     # Subtract the portion of the sphere's area that has already been covered by previous collisions
        #     portion_area -= total_portion_area

        #     # If the portion area is negative, set it to zero to avoid double-counting
        #     if portion_area < 0:
        #         portion_area = 0

        #     # Add the portion area to the total portion area
        #     total_portion_area += portion_area

        #     # Print the result for this collision
        #     print("Collision {}: Theta = {} deg, Phi = {} deg, covers {:.2f}% of the sphere's area.".format(i+1, theta_collisions[i], phi_collisions[i], portion_area))




        # # Stack the Cartesian coordinates into a matrix
        # points = np.column_stack((x, y, z))
        # # print(points, points.shape[0])

        # # Compute the Spherical Voronoi diagram
        # if points.shape[0] > 4:
        #     sv = SphericalVoronoi(points, radius=1)
        # else:
        #     pass

        # # Compute the total solid angle of the diagram
        # total_solid_angle = sv.calculate_areas().sum()
        # print(total_solid_angle)

        # # Check if the solid angle covers at least half of the sphere
        # print(ijk)
        # if coverage >= 50:
        #     print("The points cover at least half of the sphere.", coverage)
        #     bs.write("ATOM"+" "*(7-len(str(n)))+str(n)+"  O   HOH X"+" "*(4-len(str(n)))+str(n)+" "*(8-len(str(str(ijk[0]).split('.')[0])))+str(format(ijk[0], ".3f"))+" "*(4-len(str(str(ijk[1]).split('.')[0])))+str(format(ijk[1], ".3f"))+" "*(4-len(str(str(ijk[2]).split('.')[0])))+str(format(ijk[2], ".3f"))+"  1.00 30.00           O  \n")
        #     n += 1
        # else:
        #     print("The points do not cover at least half of the sphere.")

        # for start_angle in range(0, 181, 30):
        #     end_angle = start_angle + 181
        #     range_list = list(range(start_angle, end_angle, 30))
        #     if all(elem in theta_collisions for elem in range_list):
        #         for start_angle2 in range(0, 121, 30):
        #             end_angle2 = start_angle2 + 241
        #             range_list2 = list(range(start_angle2, end_angle2, 30))
        #             if all(elem in phi_collisions for elem in range_list2):
        #                 ijk = np.round(ijk, decimals=3)
        #                 print(ijk)
        #                 print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)

        #                break
                        #print('yes')
                        
    # bs.write("END                                                                             \n")    




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

sasa = ShrakeRupley2()
sasa_values = np.array(sasa.compute(coords, radii))

exposed_coords = []
exposed_radii = []
for s, c, r, sas in zip(sasa_values, coords, radii, sasa_values):
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
# print(angles_combination)
get_angles_and_atoms_3(surface_grid, angles_combination, radii, tree)


print(f"finished after {round(time() - start, 2)} seconds")


