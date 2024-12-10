from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import modules.rcparams
from numpy.linalg import norm
from numpy import cross, dot

################################# LINEAR FIT FUNCTION ######################################


def f(x, m, c):
    return m * x + c

################################# ENTER PARAMETERS ########################################


# ----------------- DATA FILE PARAMETERS -----------------
cell_info_file = 'info/box_info.txt'
filament_info_file = 'info/filament_info.txt'

# ----------------- SAMPLING PARAMETERS -----------------
sample_window_fraction = 0.02

# ----------------- EXPECTED DIFFUSION COEFFICIENT -----------------
# In units of (microns)^2 / sec
D_expected = 5

################################# READ DATA ###############################################

# Read filament info
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

try:
    data_file = 'mon_pos/mon_pos.{}.txt'.format(run_i)
    raw_data = np.loadtxt(data_file, unpack=True)
except FileNotFoundError:
    raise FileNotFoundError(
        'File index {} not found in mon_pos directory.'.format(run_i))

# --------------------------------------------------------------------------------------------

# The first row of the data file contains the time steps
t_list = raw_data[0]
num_iterations = len(t_list)

# The first row of the data file contains the time steps
t_list = raw_data[0]
num_iterations = len(t_list)

################################# ANALYZE DATA ############################################

mon_pos = np.zeros((num_iterations, num_monomers, 3))

# --------------------------------------------------------------------------------------------

# The data file contains the positions of all monomers at each time step
# Save a matrix of monomer positions for each time step
for t_i, t in enumerate(t_list):
    for m_i in range(num_monomers):
        px = raw_data[1 + 3 * m_i][t_i]
        py = raw_data[2 + 3 * m_i][t_i]
        pz = raw_data[3 + 3 * m_i][t_i]
        mon_pos[t_i, m_i] = [px, py, pz]

# Array to store radius of curvature
r_curvature = np.zeros(num_monomers - 2)
centers_of_curvature = np.zeros((num_monomers - 2, 3))

# --------------------------------------------------------------------------------------------

# Store the monomer locations at final time step

mon_pos_final = mon_pos[-1, :]

# --------------------------------------------------------------------------------------------

# Solution obtained from https://en.wikipedia.org/wiki/Circumcircle

for m_i in range(0, num_monomers - 2):
    # Calculate the vectors between the monomers
    p1 = mon_pos_final[m_i]
    p2 = mon_pos_final[m_i + 1]
    p3 = mon_pos_final[m_i + 2]

    a = p1 - p3
    b = p2 - p3

    r = norm(a) * norm(b) * norm(a - b) / (2 * norm(np.cross(a, b)))

    p_center = p3 + cross(((norm(a))**2 * b - (norm(b))**2 * a),
                          cross(a, b)) / (2 * (norm(cross(a, b)))**2)

    r_curvature[m_i] = r
    
    centers_of_curvature[m_i] = p_center

s_list = np.arange(0, num_monomers - 2) * l

plt.plot(s_list, r_curvature, label='Data', marker='o', color='black')

plt.xlabel(r'$s\,\mathrm{(nm)}$')
plt.ylabel(r'$r\,\mathrm{(nm)}$')

plt.savefig('plots/curvature.{}.pdf'.format(run_i))