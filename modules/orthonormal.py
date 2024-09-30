import numpy as np
from numpy.linalg import norm
from numpy import cross, dot


def face_vector(h):
    h = np.array(h)

    h_x, h_y, h_z = h[0], h[1], h[2]

    # if hx != 0:
    #     fx = -(1/hx) * (hy + hz)
    #     fy = 1
    #     fz = 1
    # elif hx == 0 and hy != 0:
    #     fx = 1
    #     fy = - hz / hy
    #     fz = 1
    # elif hx == 0 and hy == 0 and hz != 0:
    #     fx = 1
    #     fy = 0
    #     fz = 0

    h_x, h_y, h_z = h

    if h_x != 0 and h_y == 0 and h_z == 0:
        # h is along the x-axis
        f = (0, 1, 0)
        g = (0, 0, 1)
    elif h_x == 0 and h_y != 0 and h_z == 0:
        # h is along the y-axis
        f = (0, 0, 1)
        g = (1, 0, 0)
    elif h_x == 0 and h_y == 0 and h_z != 0:
        # h is along the z-axis
        f = (1, 0, 0)
        g = (0, 1, 0)
    else:
        # General case
        if h_x != 0 or h_y != 0 or h_z != 0:
            e_x = (1, 0, 0) if not (h_x == 1 and h_y ==
                                    0 and h_z == 0) else (0, 1, 0)
            f = cross(h, e_x)
            f = f / norm(f)
            g = cross(h, f)
            g = g / norm(g)

    # return np.array([fx, fy, fz])
    return f


def orthonormal_pair(h):
    h = np.array(h)

    f = face_vector(h)
    f = f / np.linalg.norm(f)

    g = np.cross(h, f)
    g = g / np.linalg.norm(g)

    return f, g
