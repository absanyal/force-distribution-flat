import numpy as np
from analysis_modules.cylindermath import cylinder


def make_box(box_info_filename: str):
    xlo, xhi, ylo, yhi, zlo, zhi, x_len, y_len, z_len = np.loadtxt(
        box_info_filename)

    if x_len > y_len and x_len > z_len:
        # Longest dimension is x
        radius = y_len / 2
        rA = [0, (yhi + ylo) / 2, (zhi + zlo) / 2]
        rB = [xhi, (yhi + ylo) / 2, (zhi + zlo) / 2]
    elif y_len > x_len and y_len > z_len:
        # Longest dimension is y
        radius = x_len / 2
        rA = [(xhi + xlo) / 2, 0, (zhi + zlo) / 2]
        rB = [(xhi + xlo) / 2, yhi, (zhi + zlo) / 2]
    elif z_len > x_len and z_len > y_len:
        # Longest dimension is z
        radius = x_len / 2
        rA = [(xhi + xlo) / 2, (yhi + ylo) / 2, 0]
        rB = [(xhi + xlo) / 2, (yhi + ylo) / 2, zhi]
    else:
        raise ValueError("Provided box dimensions are invalid.")

    cyl = cylinder(radius, rA, rB)
    return cyl
