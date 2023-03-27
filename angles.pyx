from scipy.spatial import cKDTree
import numpy as np
import math

def spherical_to_cartesian(r, theta, phi):
    theta = theta/180*math.pi
    phi = phi/180*math.pi
    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)
    x = round(x, 3)
    y = round(y, 3)
    z = round(z, 3)
    return [x, y, z]

def get_angles_and_atoms_3(outside_points, angles_combination, radii, tree_coords):
    #radii_inside = radii[tree_coords.query(outside_points)[1]]
    mask = tree_coords.query(outside_points)[0] <= radii[tree_coords.query(outside_points)[1]]
    
    for ijk, m in zip(outside_points, mask):
        # print(ijk)
        if m:
            print('Inside', ijk)
            continue
        theta_collisions = []
        phi_collisions = []
        for angle in angles_combination:
            t = angle[0]
            p = angle[1]
            points = np.array([np.add(np.array(ijk), np.array(spherical_to_cartesian(d, t, p))) for d in range(0, 25)])
            collisions = tree_coords.query(points)[0] <= radii[tree_coords.query(ijk)[1]]
            if np.any(collisions):
                theta_collisions.append(t)
                phi_collisions.append(p)
        #print('Theta collision:', theta_collisions, 'Phi collision:', phi_collisions)