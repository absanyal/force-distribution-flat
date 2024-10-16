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
plot_num_attached = 1

# ----------------- SMOOTHING PARAMETERS -----------------
smoothing_window = 11
detection_window = 21
threshold = 0.2
hitting_distance = 2.0
fraction_to_skip_before_recording = 0.1


################################# READ DATA ###############################################

cell = make_box(cell_info_file)
try:
    R, d, a, a1, a2, l, s1, s2, aF, aL, theta1, theta2, gamma, phi1, phi2, phi3, phi4, num_monomers, num_layers, num_total_particles, num_linkers, num_bonds, num_angles = np.loadtxt(
        filament_info_file)
except FileNotFoundError:
    raise FileNotFoundError(
        'Please provide the correct path to the filament info file.')

num_monomers = int(num_monomers)
num_linkers = int(num_linkers)

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

t_list = raw_data[0]
num_iterations = len(t_list)
recording_start_index = int(fraction_to_skip_before_recording * num_iterations)

################################# ANALYSIS ###############################################

# ----------------- DISTANCE FROM SURFACE -----------------

distances = np.zeros((num_iterations, num_linkers))
distances_smooth = np.zeros_like(distances)
detection_smooth = np.zeros_like(distances)

# Save the distances from the surface for each linker at each time step
for t_i, t in enumerate(t_list):
    for l_i in range(num_linkers):
        x = raw_data[1 + 3*l_i][t_i]
        y = raw_data[2 + 3*l_i][t_i]
        z = raw_data[3 + 3*l_i][t_i]

        rP = [x, y, z]
        dist = distance_from_surface(cell, rP)
        distances[t_i, l_i] = dist

# Smoothen the distances
for l_i in range(num_linkers):
    distances_smooth[:, l_i] = moving_average(
        distances[:, l_i], smoothing_window, padding='edge')
    detection_smooth[:, l_i] = moving_average(
        distances[:, l_i], detection_window, padding='edge')


# Plot the distances
if plot_traces:

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    for l_i in range(num_linkers):
        ax.plot(t_list, distances_smooth[:, l_i], label='Linker {}'.format(
            l_i+1), linewidth=0.5)

    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'$\Delta s_{\mathrm{linker}}$ (nm)')
    ax.set_xlim([t_list[0], t_list[-1]])

    ax.axhline(hitting_distance, color='black',
               linestyle='--', label='Hitting distance')

    plt.savefig('plots/linker_distances.{}.pdf'.format(run_i), dpi=300)

# ----------------- Number of linkers attached -----------------

num_attached = np.zeros(num_iterations)

for t_i in range(num_iterations):
    for l_i in range(num_linkers):
        s = abs(detection_smooth[t_i, l_i] -
                hitting_distance) / hitting_distance
        if s < threshold:
            num_attached[t_i] += 1

avg_num_attached = np.mean(num_attached[recording_start_index:])

fraction_attached = num_attached / num_linkers
avg_fraction_attached = avg_num_attached / num_linkers

if plot_num_attached:

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    # ax.plot(t_list, num_attached, label='Number of attached linkers', linewidth=1.0, color='k')
    ax.plot(t_list, fraction_attached, label='Fraction of attached linkers', linewidth=0.5, color='k')

    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'$N_{\mathrm{attached}} / N_{\mathrm{linkers}}$')
    ax.set_xlim([t_list[0], t_list[-1]])

    # ax.axhline(avg_num_attached, color='r', linestyle='--',
    #            linewidth=0.5, label='Average: {:.2f}'.format(avg_num_attached))
    ax.axhline(avg_fraction_attached, color='r', linestyle='--', linewidth=0.5, label='Average fraction: {:.2f}'.format(avg_fraction_attached))

    ax.legend()
    plt.savefig('plots/num_attached.{}.pdf'.format(run_i), dpi=300)
    
    with open('data/num_attached.{}.txt'.format(run_i), 'w') as f:
        f.write('# avg_num_attached \t avg_fraction_attached\n')
        f.write('{} \t {}'.format(avg_num_attached, avg_fraction_attached))
