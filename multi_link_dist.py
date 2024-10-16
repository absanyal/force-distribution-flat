from sys import argv
import numpy as np
from analysis_modules.cylindermath import cylinder, moving_average, distance_from_surface
from analysis_modules.create_cell import make_box
import matplotlib.pyplot as plt
import modules.rcparams

cell_info_file = 'info/box_info.txt'
filament_info_file = 'info/filament_info.txt'

smoothing_window = 51

cell = make_box(cell_info_file)

R, d, a, a1, a2, l, s1, s2, aF, aL, theta1, theta2, gamma, phi1, phi2, phi3, phi4, num_monomers, num_layers, num_total_particles, num_linkers, num_bonds, num_angles = np.loadtxt(
    filament_info_file)
num_monomers = int(num_monomers)
num_linkers = int(num_linkers)

if len(argv) < 1:
    raise ValueError('Please provide a run index.')
elif len(argv) > 2:
    raise ValueError('Too many arguments provided.')

run_i = int(argv[1])
if run_i < 0:
    raise ValueError('Run index must be a non-negative integer.')

data_file = 'link_pos/link_pos.{}.txt'.format(run_i)

raw_data = np.loadtxt(data_file, unpack=True)

t_list = raw_data[0]
num_iterations = len(t_list)

distances = np.zeros((num_iterations, num_linkers))
distances_smooth = np.zeros_like(distances)

for t_i, t in enumerate(t_list):
    for l_i in range(num_linkers):
        x = raw_data[1 + 3*l_i][t_i]
        y = raw_data[2 + 3*l_i][t_i]
        z = raw_data[3 + 3*l_i][t_i]
        
        rP = [x, y, z]
        dist = distance_from_surface(cell, rP)
        distances[t_i, l_i] = dist

for l_i in range(num_linkers):
    distances_smooth[:, l_i] = moving_average(distances[:, l_i], smoothing_window, padding='edge')
        

for l_i in range(num_linkers):
    plt.plot(t_list, distances_smooth[:, l_i], label='Linker {}'.format(l_i+1), linewidth=0.5)

plt.xlabel('Time')
plt.ylabel('Distance from surface')
plt.xlim([t_list[0], t_list[-1]])
plt.show()