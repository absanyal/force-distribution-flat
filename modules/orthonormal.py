import numpy as np

def face_vector(h):
    h = np.array(h)
    
    hx, hy, hz = h[0], h[1], h[2]
    
    if hx != 0:
        fx = -(1/hx) * (hy + hz)
        fy = 1
        fz = 1
    elif hx == 0 and hy != 0:
        fx = 1
        fy = - hz / hy
        fz = 1
    elif hx == 0 and hy == 0 and hz != 0:
        fx = 1
        fy = 0
        fz = 0
    
    return np.array([fx, fy, fz])

def orthonormal_pair(h):
    h = np.array(h)
    
    f = face_vector(h)
    f = f / np.linalg.norm(f)
    
    g = np.cross(h, f)
    g = g / np.linalg.norm(g)
    
    return f, g