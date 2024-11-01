from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import modules.rcparams

################################# ENTER PARAMETERS ########################################

# ----------------- DATA FILE PARAMETERS -----------------
cell_info_file = 'info/box_info.txt'
filament_info_file = 'info/filament_info.txt'

# ----------------- PLOT TOGGLES -------------------------

plot_e2e_distance = True

# ----------------- SMOOTHING PARAMETERS -----------------
fraction_to_skip_before_recording = 0.05

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
lc = num_monomers * a # Contour length of the filament

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

# Averages are calculated after this many iterations have passed
recording_start_index = int(fraction_to_skip_before_recording * num_iterations)

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

# --------------------------------------------------------------------------------------------

# Calculate the end-to-end distance
e2e_dist = np.zeros(num_iterations)

for t_i in range(num_iterations):
    mon_start = mon_pos[t_i, 0]
    mon_end = mon_pos[t_i, -1]
    
    e2e_dist[t_i] = np.linalg.norm(mon_end - mon_start)

avg_e2e_dist = np.mean(e2e_dist[recording_start_index:])
    
if plot_e2e_distance:
    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))
    
    ax.plot(t_list, e2e_dist, color='black', lw=1, label='End-to-end distance')
    ax.axhline(avg_e2e_dist, color='red', lw=1, ls='--', label='Average: {:.2f}'.format(avg_e2e_dist))
    
    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'$R_{\mathrm{e-e}}$')
    
    ax.legend()
    
    plt.savefig('plots/e2e_distance.{}.pdf'.format(run_i), dpi=300)


# --------------------------------------------------------------------------------------------

# Calculate the tangent vectors
tangent_vectors_list = np.zeros((num_iterations, num_monomers - 1, 3))

for t_i in range(num_iterations):
    for m_i in range(num_monomers - 1):
        tangent_vectors_list[t_i, m_i] = mon_pos[t_i, m_i + 1] - mon_pos[t_i, m_i]