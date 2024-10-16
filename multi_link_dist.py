from sys import argv
import numpy as np
from analysis_modules.cylindermath import cylinder, moving_average, distance_from_surface
from analysis_modules.create_cell import make_box

cell_info_file = 'info/box_info.txt'
filament_info_file = 'info/filament_info.txt'

cell = make_box(cell_info_file)
R, d, a, a1, a2, l, s1, s2, aF, aL, theta1, theta2, gamma, phi1, phi2, phi3, phi4, num_monomers, num_layers, num_total_particles, num_linkers, num_bonds, num_angles = np.loadtxt(
    filament_info_file)

if len(argv) < 1:
    raise ValueError('Please provide a run index.')
elif len(argv) > 2:
    raise ValueError('Too many arguments provided.')

run_index = int(argv[1])

if run_index < 0:
    raise ValueError('Run index must be a non-negative integer.')
