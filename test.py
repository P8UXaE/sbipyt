import numpy as np
import itertools
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
    # return [x, y, z] # This if used with @njit
    return np.array([x, y, z]) # This if not used with @njit

theta = list(range(0, 361, 45))
phi = list(range(0, 361, 45))
angles_combination = np.array(list(itertools.product(theta, phi)))
print(angles_combination)
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
print(unique_angles)