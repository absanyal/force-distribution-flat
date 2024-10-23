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

# ----------------- SMOOTHING PARAMETERS -----------------
smoothing_window = 11
detection_window = 51
threshold = 0.2
hitting_distance = 2.0
fraction_to_skip_before_recording = 0.05

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

# --------------------------------------------------------------------------------------------

# The first row of the data file contains the time steps
t_list = raw_data[0]
num_iterations = len(t_list)

# Averages are calculated after this many iterations have passed
recording_start_index = int(fraction_to_skip_before_recording * num_iterations)

################################# ANALYSIS ###############################################

# ----------------- Assign segments to linkers -----------------

segment_id = []
linker_distribution = np.loadtxt('info/linker_distribution.txt')
for i in range(num_monomers):
    if linker_distribution[i] == 1:
        segment_id.append(i)


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
            l_i+1), linewidth=0.5)

    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'$\Delta s_{\mathrm{linker}}$ (nm)')
    ax.set_xlim([t_list[0], t_list[-1]])
    ax.set_ylim(bottom=0)

    # ax.axhline(hitting_distance, color='black',
    #            linestyle='--', label='Hitting distance')

    plt.savefig('plots/linker_distances.{}.pdf'.format(run_i), dpi=300)
    
    plt.clf()
    plt.cla()
    ax.cla()

# ----------------- PROXIMITY TO SURFACE ----------------------

target_list = distances_list

proximity_list = np.zeros((num_iterations, num_monomers))
minimum_nonzero_proximity = abs(
    cell.radius - np.max(distances_list[np.nonzero(target_list)])) / cell.radius
maximum_proximity = abs(cell.radius - np.min(target_list)) / cell.radius

        

# Calculate the proximity of each linker to the surface
for t_i in range(num_iterations):
    for l_i in range(num_linkers):
        s_i = segment_id[l_i]
        proximity = abs(cell.radius - target_list[t_i, l_i]) / cell.radius
        proximity_list[t_i, s_i] = proximity

if plot_proximity:

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    im = ax.imshow(proximity_list.T, aspect='auto', cmap='plasma', origin='lower', interpolation='None', extent=[
              t_list[0], t_list[-1], 1, num_monomers+1], vmin=minimum_nonzero_proximity, vmax=maximum_proximity)

    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'Proximity to surface')

    plt.savefig('plots/proximity.{}.pdf'.format(run_i), dpi=300)
    
    plt.clf()
    plt.cla()
    ax.cla()

# ----------------- ATTACHMENT STATUS PER LINKER ----------------------

target_list = distances_smooth_list

attachment_status = np.zeros((num_iterations, num_monomers))
for t_i in range(num_iterations):
    for l_i in range(num_linkers):
        s_i = segment_id[l_i]
        s = abs(target_list[t_i, l_i] -
                hitting_distance) / hitting_distance
        if s <= threshold:
            attachment_status[t_i, s_i] = 1.0

if plot_attachment_status:
    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))
    
    ax.imshow(attachment_status.T, aspect='auto', cmap='binary', origin='lower', extent=
                [t_list[0], t_list[-1], 1, num_monomers+1])
    
    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'Attachment status')
    
    plt.savefig('plots/attachment_status.{}.pdf'.format(run_i), dpi=300)
    
    plt.clf()
    plt.cla()
    ax.cla()

# ----------------- ATTACHED SEGMENT LENGTH ----------------------

attached_segment_length = np.zeros(num_iterations)

for t_i in range(num_iterations):
    attached_segement_indices = np.where(attachment_status[t_i] == 1.0)[0]
    if len(attached_segement_indices) != 0:
        first_index = attached_segement_indices[0]
        last_index = attached_segement_indices[-1]
        attached_segment_length[t_i] = last_index - first_index + 1
    else:
        attached_segment_length[t_i] = 0

average_attached_segment_length = np.mean(attached_segment_length[recording_start_index:])
maximal_attached_segment_length = np.max(attached_segment_length)
    
if plot_attached_segment_length:
    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))
    
    ax.plot(t_list, attached_segment_length, linewidth=1.0, color='k')
    
    ax.axhline(average_attached_segment_length, color='r', linestyle='--', linewidth=0.5, label='Average: {:.2f}'.format(average_attached_segment_length))
    ax.axhline(maximal_attached_segment_length, color='b', linestyle='--', linewidth=0.5, label='Maximum: {}'.format(maximal_attached_segment_length))
    
    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'Attached segment length')
    ax.set_xlim([t_list[0], t_list[-1]])
    
    ax.legend(loc='upper right')
    
    plt.savefig('plots/attached_segment_length.{}.pdf'.format(run_i), dpi=300)
    
    plt.clf()
    plt.cla()
    ax.cla()

# ----------------- NUMBER OF ATTACHED LINKERS ----------------------

target_list = detection_smooth_list

# Count the number of linkers attached at each time step
num_attached = np.zeros(num_iterations)

for t_i in range(num_iterations):
    for l_i in range(num_linkers):
        s = abs(target_list[t_i, l_i] -
                hitting_distance) / hitting_distance
        if s < threshold:
            num_attached[t_i] += 1

# Calculate the average number of attached linkers
avg_num_attached = np.mean(num_attached[recording_start_index:])

fraction_attached = num_attached / num_linkers
avg_fraction_attached = avg_num_attached / num_linkers

# Maximum number of attached linkers
max_num_attached = np.max(num_attached)

# Plot the number of attached linkers
if plot_num_attached:

    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    ax.plot(t_list, num_attached, linewidth=1.0, color='k')
    # ax.plot(t_list, fraction_attached, label='Fraction of attached linkers', linewidth=0.5, color='k')

    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'$N_{\mathrm{attached}}$')
    ax.set_xlim([t_list[0], t_list[-1]])

    ax.axhline(avg_num_attached, color='r', linestyle='--',
               linewidth=0.5, label='Average: {:.2f}'.format(avg_num_attached))
    ax.axhline(max_num_attached, color='b', linestyle='--',
               linewidth=0.5, label='Maximum: {}'.format(max_num_attached))

    # ax.axhline(avg_fraction_attached, color='r', linestyle='--', linewidth=0.5, label='Average fraction: {:.2f}'.format(avg_fraction_attached))

    # ax.grid(axis='y', linestyle='--', linewidth=0.5, color='gray', which='major')

    ax.legend(loc='upper right')
    plt.savefig('plots/num_attached.{}.pdf'.format(run_i), dpi=300)
    
    plt.clf()
    plt.cla()
    ax.cla()

    with open('data/num_attached.{}.txt'.format(run_i), 'w') as f:
        f.write('# avg_num_attached \t avg_fraction_attached\n')
        f.write('{} \t {}'.format(avg_num_attached, avg_fraction_attached))
