from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import modules.rcparams
from scipy.optimize import curve_fit
from numpy import cos, sin, log

################################# ENTER PARAMETERS ########################################

# ----------------- DATA FILE PARAMETERS -----------------
cell_info_file = 'info/box_info.txt'
filament_info_file = 'info/filament_info.txt'

# ----------------- PLOT TOGGLES -------------------------

plot_e2e_distance = 1
plot_e2e_hist = 1
plot_correlations = 1
plot_correlations_linearized = 1

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

# Calculate the end-to-end distance
e2e_dist = np.zeros(num_iterations)

for t_i in range(num_iterations):
    mon_start = mon_pos[t_i, 0]
    mon_end = mon_pos[t_i, -1]

    e2e_dist[t_i] = np.linalg.norm(mon_end - mon_start)

avg_e2e_dist = np.mean(e2e_dist[recording_start_index:])
expected_e2e_dist = 2 * R * sin(lc / (2 * R))
spread_e2e_dist = np.std(e2e_dist[recording_start_index:])
lower_bound = avg_e2e_dist - spread_e2e_dist
upper_bound = avg_e2e_dist + spread_e2e_dist

if plot_e2e_distance:
    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    ax.plot(t_list, e2e_dist, color='black', lw=1, label='End-to-end distance')
    ax.axhline(avg_e2e_dist, color='red', lw=3, ls='--',
               label=r'Average: ${:.2f} \pm {:.2f}$ nm'.format(avg_e2e_dist, spread_e2e_dist))
    ax.axhline(expected_e2e_dist, color='g', lw=3, ls='--',
               label='Expected: {:.2f} nm'.format(expected_e2e_dist))
    ax.fill_between(t_list, lower_bound, upper_bound, color='red', alpha=0.2)

    ax.set_xlabel(r'$t/\tau$')
    ax.set_ylabel(r'$R_{\mathrm{e-e}}$')

    ax.set_xlim(t_list[0], t_list[-1])

    ax.legend()

    plt.savefig('plots/e2e_distance.{}.pdf'.format(run_i), dpi=300)

if plot_e2e_hist:
    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    ax.hist(e2e_dist[recording_start_index:], bins='auto',
            color='blue', edgecolor='none', rwidth=0.8, density=True)
    ax.axvline(avg_e2e_dist, color='red', lw=3, ls='--',
               label=r'Average: ${:.2f} \pm {:.2f}$ nm'.format(avg_e2e_dist, spread_e2e_dist))
    ax.axvline(expected_e2e_dist, color='g', lw=3, ls='--',
               label='Expected: {:.2f} nm'.format(expected_e2e_dist))

    ax.axvline(lower_bound, color='red', lw=0.5, ls='--')
    ax.axvline(upper_bound, color='red', lw=0.5, ls='--')

    ax.set_xlabel(r'$R_{\mathrm{e-e}}$')
    ax.set_ylabel('Frequency')

    ax.legend()

    plt.savefig('plots/e2e_hist.{}.pdf'.format(run_i), dpi=300)

# --------------------------------------------------------------------------------------------

# Calculate the tangent vectors
tangent_vectors_list = np.zeros((num_iterations, num_monomers - 1, 3))
correlations = np.zeros((num_iterations, num_monomers - 1))
average_correlations = np.zeros(num_monomers - 1)
# correlation_errors = np.zeros(num_monomers - 1)

average_correlations_linearized = np.zeros_like(average_correlations)

s_list = np.zeros_like(average_correlations)
for s in range(num_monomers - 1):
    s_list[s] = s * a

for t_i in range(num_iterations):
    for m_i in range(num_monomers - 1):
        tangent_vectors_list[t_i, m_i] = mon_pos[t_i,
                                                 m_i + 1] - mon_pos[t_i, m_i]

for t_i in range(num_iterations):
    for vi in range(len(tangent_vectors_list[t_i])):
        tangent_0 = tangent_vectors_list[t_i, 0]
        tangent_s = tangent_vectors_list[t_i, vi]
        correlations[t_i, vi] = np.dot(
            tangent_0, tangent_s) / (np.linalg.norm(tangent_0) * np.linalg.norm(tangent_s))

for m_i in range(num_monomers - 1):
    average_correlations[m_i] = np.mean(correlations[:, m_i])
    # correlation_errors[m_i] = np.std(correlations[:, m_i])

# for m_i in range(num_monomers - 1):
#     average_correlations[m_i] = average_correlations[m_i] / ( cos((s_list[m_i]) / R) )

for m_i in range(num_monomers - 1):
    s = s_list[m_i]
    average_correlations_linearized[m_i] = log(average_correlations[m_i] / cos(s / R))

# --------------------------------------------------------------------------------------------

# Fitting the data


def correlation_function(s, l_p):
    return np.exp(-(s) / l_p) * cos((s) / R)

# def correlation_function(s, l_p, alpha):
#     return np.exp(-(s) / l_p) * cos(s/alpha)

def linear(s, lp):
    return -s / lp


popt, pcov = curve_fit(correlation_function, s_list, average_correlations)


# lp_fit, alpha_fit = popt
# err_lp, err_alpha = np.sqrt(np.diag(pcov))

lp_fit, = popt
err_lp, = np.sqrt(np.diag(pcov))

# print('Persistence length: {:.2f} +/- {:.2f} nm'.format(lp_fit, err_lp))
# print('Alpha: {:.2f} +/- {:.2f}'.format(alpha_fit, err_alpha))

fitting_s = np.linspace(s_list[0], s_list[-1], 100)
fitting_correlations = correlation_function(s_list, *popt)

# --------------------------------------------------------------------------------------------

if plot_correlations:
    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    ax.plot(s_list, average_correlations, color='black', lw=0,
            label='Simulation', marker='o', markersize=3)

    # ax.plot(s_list, fitting_correlations, color='red', marker='o', markersize = 2, lw=1, ls='--', label=r'$l_p = {:.2f} \pm {:.2f}\,\mathrm{{nm}}\\\,\alpha = {:.2f} \pm {:.2f}\,\mathrm{{nm}}$'.format(lp_fit, err_lp, alpha_fit, err_alpha))

    ax.plot(s_list, fitting_correlations, color='red', marker='o', markersize=2, lw=1,
            ls='--', label=r'$l_p = {:.2f} \pm {:.2f}\,\mathrm{{nm}}$'.format(lp_fit, err_lp))

    # ax.errorbar(s_list, average_correlations, yerr=correlation_errors,
    #             fmt='none', ecolor='blue', capsize=2, capthick=0.5, elinewidth=0.5)

    ax.set_xlabel(r'$s\,[\mathrm{nm}]$')
    ax.set_ylabel(r'$\langle \hat{t}_0 \cdot \hat{t}_{s} \rangle$')

    # ax.set_xscale('log')

    # plt.title(r'$\exp(-s/l_p) \cos(s / R_0)$ at $R_0 = {:.2f}\,\mathrm{{nm}}$'.format(R))
    plt.title(r'$R_0 = {:.2f}\,\mathrm{{nm}}$'.format(R))

    plt.legend()

    plt.savefig('plots/correlations.{}.pdf'.format(run_i), dpi=300)

if plot_correlations_linearized:
    
    popt, pcov = curve_fit(linear, s_list, average_correlations_linearized)
    
    lp_fit, = popt
    err_lp, = np.sqrt(np.diag(pcov))
    
    fitting_correlations_linearized = linear(s_list, *popt)
    
    fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

    ax.plot(s_list, average_correlations_linearized, color='black', lw=0,
            label='Simulation', marker='o', markersize=3)
    
    ax.plot(s_list, fitting_correlations_linearized, color='red', marker='o', markersize=2, lw=1,
            ls='--', label=r'$l_p = {:.2f} \pm {:.2f}\,\mathrm{{nm}}$'.format(lp_fit, err_lp))

    ax.set_xlabel(r'$s\,[\mathrm{nm}]$')
    ax.set_ylabel(r'$\ln\left(\langle \hat{t}_0 \cdot \hat{t}_{s} \rangle / \cos(s/R_0)\right)$')

    plt.title(r'$R_0 = {:.2f}\,\mathrm{{nm}}$'.format(R))

    plt.legend()

    plt.savefig('plots/correlations_linearized.{}.pdf'.format(run_i), dpi=300)
