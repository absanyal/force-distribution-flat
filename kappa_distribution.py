from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import modules.rcparams
from scipy.optimize import curve_fit
from numpy import cos, sin, log, dot
from numpy.linalg import norm
import scipy.stats as stats

################################# ENTER PARAMETERS ########################################

# ----------------- DATA FILE PARAMETERS -----------------
cell_info_file = 'info/box_info.txt'
filament_info_file = 'info/filament_info.txt'

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
lc = a * (num_monomers - 1)

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

# Calculate the tangent vectors
tangent_vectors_list = np.zeros((num_iterations, num_monomers - 1, 3))
dtheta_ds_list = []

s_list = np.zeros(num_monomers - 1)
for s in range(num_monomers - 1):
    s_list[s] = s * a

for t_i in range(num_iterations):
    for m_i in range(num_monomers - 1):
        tangent_vectors_list[t_i, m_i] = mon_pos[t_i,
                                                 m_i + 1] - mon_pos[t_i, m_i]

for t_i in range(num_iterations):
    for v_i in range(num_monomers - 2):
        for v_j in range(num_monomers - 2):
            if (v_j - v_i) == 1:
                tangent1 = tangent_vectors_list[t_i, v_i]
                tangent2 = tangent_vectors_list[t_i, v_j]
                
                # dot_product = dot(tangent1, tangent2) / (norm(tangent1) * norm(tangent2))
                # dtheta = np.arccos(dot_product)
                
                theta1 = np.arctan2(tangent1[1], tangent1[0])
                theta2 = np.arctan2(tangent2[1], tangent2[0])
                dtheta = theta2 - theta1
                
                s1 = s_list[v_i]
                s2 = s_list[v_j]
                
                if s2 - s1 != 0:
                    dtheta_ds = dtheta / (s2 - s1)
                    dtheta_ds_list.append(dtheta_ds)

# for t_i in range(num_iterations):
#     for v_i in range(num_monomers - 2):
#         tangent1 = tangent_vectors_list[t_i, v_i]
#         tangent2 = tangent_vectors_list[t_i, v_i + 1]
#         dot_product = dot(tangent1, tangent2) / (norm(tangent1) * norm(tangent2))
#         dtheta = np.arccos(dot_product)
        
#         s1 = s_list[v_i]
#         s2 = s_list[v_i + 1]
        
#         if s2 - s1 == 0:
#             continue
#         else:
#             dtheta_ds = dtheta / (s2 - s1)
#             dtheta_ds_list.append(dtheta_ds)

# --------------------------------------------------------------------------------------------

# Calculate mean curvature

# mean_curvature = np.mean(dtheta_ds_list)
# mean_radius = 1 / mean_curvature

# print('Mean radius of curvature: {:.4f} nm'.format(mean_radius))

# n, bins, patches = plt.hist(dtheta_ds_list, bins='auto', color='blue', alpha=0.7, edgecolor='none', rwidth=0.8, density=True)

# elem = np.argmax(n)

# most_freq_kappa = (bins[elem] + bins[elem + 1]) / 2
# most_freq_radius = 1 / most_freq_kappa

# print('Most frequent radius: {:.4f} nm'.format(most_freq_radius))

ds = abs(s_list[1] - s_list[0])

sigma = np.std(dtheta_ds_list)
lp = ((1 / sigma) ** 2) * (1/(ds))

print('Standard deviation: {:.4f}'.format(sigma))
print('Persistence length: {:.4f} nm'.format(lp))


# --------------------------------------------------------------------------------------------

plt.figure(constrained_layout=True)

plt.hist(dtheta_ds_list, bins='auto', color='blue', alpha=0.5, edgecolor='none', rwidth=0.85, density=True)

# plt.axvline(mean_curvature, color='red', lw=2, ls='--', label=r'Mean radius: ${:.2f}$ nm'.format(mean_radius))

# plt.axvline(most_freq_kappa, color='g', lw=2, ls='--', label=r'Most frequent radius: ${:.2f}$ nm'.format(most_freq_radius))

plt.xlabel(r'$\frac{\mathrm{d}\theta}{\mathrm{d} s}\,(\mathrm{rad\,{nm}^{-1}})$')
plt.ylabel('Frequency')

# plt.xlim(-0.5, 0.5)

# plt.legend()

plt.savefig('plots/kappa_distribution.pdf', format='pdf')