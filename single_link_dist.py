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

# ----------------- PLOT TOGGLES ------------------------

plot_traces = 1
plot_proximity = 1
plot_num_attached = 1
plot_attachment_status = 1
plot_attached_segment_length = 1
plot_is_attached = 1
plot_attach_detach_intervals = 1

# ----------------- SMOOTHING PARAMETERS -----------------
smoothing_window = 11
detection_window = 51
threshold = 0.2
hitting_distance = 2.0
fraction_to_skip_before_recording = 0.05
attachment_number_threshold = 2 # Number of linkers to be attached for a segment to be considered attached

if smoothing_window % 2 == 0:
    raise ValueError('Smoothing window must be an odd number.')

if detection_window % 2 == 0:
    raise ValueError('Detection window must be an odd number.')

if threshold < 0 or threshold > 1:
    raise ValueError('Threshold must be between 0 and 1.')

if hitting_distance < 0:
    raise ValueError('Hitting distance must be non-negative.')

if fraction_to_skip_before_recording < 0 or fraction_to_skip_before_recording > 1:
    raise ValueError(
        'Fraction of iterations to skip before recording must be between 0 and 1.')


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

try:
    data_file = 'link_pos/link_pos.{}.txt'.format(run_i)
    raw_data = np.loadtxt(data_file, unpack=True)
except FileNotFoundError:
    raise FileNotFoundError(
        'File index {} not found in link_pos directory.'.format(run_i))

try:
    assert(num_monomers == 1)
except AssertionError:
    raise ValueError('This script is only for single monomer simulations.')

# --------------------------------------------------------------------------------------------

# The first row of the data file contains the time steps
t_list = raw_data[0]
num_iterations = len(t_list)

# Averages are calculated after this many iterations have passed
recording_start_index = int(fraction_to_skip_before_recording * num_iterations)

# Bounds
r_max = 2.1251
r_cutoff = 3.0

# ----------------- DISTANCE FROM SURFACE -----------------

distances_list = np.zeros((num_iterations, num_linkers))
distances_smooth_list = np.zeros_like(distances_list)
detection_smooth_list = np.zeros_like(distances_list)

# Save the distances from the surface for each linker at each time step
for t_i, t in enumerate(t_list):
    for l_i in range(num_linkers):
        x = raw_data[1 + 3*l_i][t_i]
        y = raw_data[2 + 3*l_i][t_i]
        z = raw_data[3 + 3*l_i][t_i]

        rP = [x, y, z]
        dist = distance_from_surface(cell, rP)
        distances_list[t_i, l_i] = dist

# Smoothen the distances
for l_i in range(num_linkers):
    distances_smooth_list[:, l_i] = moving_average(
        distances_list[:, l_i], smoothing_window, padding='edge')
    detection_smooth_list[:, l_i] = moving_average(
        distances_list[:, l_i], detection_window, padding='edge')


# Plot the distances
if plot_traces:

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    for l_i in range(num_linkers):
        ax.plot(t_list, distances_smooth_list[:, l_i], label='Linker {}'.format(
            l_i+1), linewidth=0.5, color='black')

    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'$\Delta s_{\mathrm{linker}}$ (nm)')
    ax.set_xlim([t_list[0], t_list[-1]])
    ax.set_ylim(bottom=1.5, top=r_cutoff)

    ax.axhline(r_max, color='red', linestyle='--', label=r'$r_{{\mathrm{{max}}}} = {:.2f}\,\mathrm{{nm}}$'.format(r_max))
    
    # ax.axhline(hitting_distance, color='black',
    #            linestyle='--', label='Hitting distance')

    plt.legend()
    
    plt.savefig('plots/single_linker_distances.{}.pdf'.format(run_i), dpi=300)
    
    plt.clf()
    plt.cla()
    ax.cla()

values_counter = 0
with open('r_histogram.dat', 'w') as f:
    for s_i, s in enumerate(distances_smooth_list[:, 0]):
        if s < r_cutoff:
            f.write('{}\n'.format(s))
            values_counter += 1

print('Values written:', values_counter)