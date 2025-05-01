from sys import argv
import numpy as np
from analysis_modules.cylindermath import cylinder, moving_average, distance_from_surface
from analysis_modules.create_cell import make_box
import matplotlib.pyplot as plt
import modules.rcparams

################################# ENTER PARAMETERS ########################################

# ----------------- DATA FILE PARAMETERS -----------------
cell_info_file = 'info/box_info.txt'
filament_info_file = 'info/filament_info.txt'

# ----------------- POTENTIAL PARAMETERS -----------------
epsilon = 5000
sigma = 2.1
r_cutoff = 30

kBT0 = 310

################################# READ DATA ###############################################

# Read the box dimensions to create the cell, read filament info
cell = make_box(cell_info_file)
try:
    R, d, a, a1, a2, l, s1, s2, aF, aL, theta1, theta2, gamma, phi1, phi2, phi3, phi4, num_monomers, num_layers, num_total_particles, num_linkers, num_bonds, num_angles = np.loadtxt(
        filament_info_file)
except FileNotFoundError:
    raise FileNotFoundError(
        'Please provide the correct path to the filament info file.')

num_monomers = int(num_monomers)
num_linkers = int(num_linkers)

# --------------------------------------------------------------------------------------------
# The index of the run being analyzed is passed as an argument to the script
# Checking if the correct number of arguments are provided and if the run index is an integer
# Also checking if the file exists for the run index provided

if len(argv) == 1:
    raise ValueError('Please provide a run index.')
elif len(argv) > 2:
    args_provided = len(argv) - 1
    raise ValueError('Expected 1 argument, got {}.'.format(args_provided))

try:
    run_i = int(argv[1])
except ValueError:
    raise ValueError('Run index must be an integer.')

if run_i < 0:
    raise ValueError('Run index must be a non-negative integer.')

# try:
#     data_file = 'link_pos/link_pos_vertices.{}.txt'.format(run_i)
#     raw_data = np.loadtxt(data_file, unpack=True)
# except FileNotFoundError:
#     raise FileNotFoundError(
#         'File index {} not found in link_pos directory.'.format(run_i))

try:
    energy_file = 'thermo/energy.{}.txt'.format(run_i)
    t_e_list, temperature, E_k, E_p, E_tot = np.loadtxt(energy_file, unpack=True)
except FileNotFoundError:
    raise FileNotFoundError(
        'File index {} not found in thermo directory.'.format(run_i))

# --------------------------------------------------------------------------------------------

# The first row of the data file contains the time steps
# t_list = raw_data[0]
num_iterations = len(t_e_list)

R_cell = cell.radius

# --------------------------------------------------------------------------------------------

################################# LOAD DATA ###############################################

# Vertex data

num_vertices = num_monomers * 4

# The data file contains the positions of the vertices of the links

# The data is stored in the following format:
# t x1 y1 z1 x2 y2 z2 ... xN yN zN

distances_list = np.zeros((num_iterations, num_vertices * 3))

linker_potential_energy_list = np.zeros(num_iterations)

def LJ93(r):
    return epsilon * ( (2/15) * (sigma / r)**9 - (sigma / r)**3 ) + epsilon * ( (2/15) * (sigma / r_cutoff)**9 - (sigma / r_cutoff)**3 )

initial_total_energy = E_tot[0]
final_total_energy = E_tot[-1]

# --------------------------------------------------------------------------------------------

E_bend = final_total_energy - initial_total_energy

print("R0 = {:.2f}".format(R))
print
print("delta E_bend = {:.4f} kBT".format(E_bend / kBT0))