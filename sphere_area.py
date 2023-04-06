import numpy as np
import itertools

theta_collision= [45, 45, 45, 45, 90, 90, 90, 90, 135, 135, 135]
phi_collision= [0, 45, 90, 315, 0, 45, 90, 315, 0, 45, 315]

# theta_collision = [0, 45, 45, 45, 45, 45, 45, 45, 45, 90, 90, 90, 90, 135, 135, 135, 135, 180]
# phi_collision = [0, 0, 45, 90, 135, 180, 225, 270, 315, 0, 45, 90, 135, 45, 90, 135, 180, 0]


delta_angles = np.array([[-45, 0],[45, 0],[0, -45],[0, 45]])

pairs = np.array([(t, p) for t, p in zip(theta_collision, phi_collision)])

print(pairs)


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


# for th in pairs:
#     print('-------', th)
#     neighbors = np.array([]).reshape(0, 2) 
#     for delta in delta_angles:
#         check_coord = np.add(th, delta)
#         if check_coord[0] < -0:
#             check_coord[0] = check_coord[0]+45
#         if check_coord[0] > 180:
#             check_coord[0] = check_coord[0]-45
#         if check_coord[1] < -0:
#             check_coord[1] = check_coord[1]+360
#         if check_coord[1] >= 360:
#             check_coord[1] = check_coord[1]-360
#         if np.any(np.all(pairs == check_coord, axis=1)):
#             neighbors = np.vstack([neighbors, pairs[np.all(pairs == check_coord, axis=1)]])
#     if len(neighbors) > 1:
#         print(neighbors)

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
    print(surface)
            

print(surface)





