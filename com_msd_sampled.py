from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import modules.rcparams
from scipy.optimize import curve_fit
from numpy import cos, sin, log

from scipy.optimize import curve_fit

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

com_pos = np.zeros((num_iterations, 3))

# Calculate the center of mass of the system at each time step
for t_i in range(num_iterations):
    for m_i in range(num_monomers):
        com_pos[t_i] += mon_pos[t_i, m_i]
    com_pos[t_i] /= num_monomers

com_displacement = np.zeros((num_iterations, 3))
for t_i in range(num_iterations):
    com_displacement[t_i] = com_pos[t_i] - com_pos[0]

t_max = len(t_list)

sample_window = int(t_max * sample_window_fraction)
print("Total iterations: {}".format(t_max))
print("Sample window: {} iterations".format(sample_window))

dx = com_displacement[:, 0]
dy = com_displacement[:, 1]
dz = com_displacement[:, 2]

ds_sq_avg_sampled = np.zeros(sample_window)
t0_iter_list = np.arange(0, t_max-sample_window)
print("Samples taken: {}".format(len(t0_iter_list)))

t_shortened = t_list[:sample_window]

for t0_i in t0_iter_list:
    init_dx = dx[t0_i]
    init_dy = dy[t0_i]
    init_dz = dz[t0_i]
    
    ds_sq_avg_sampled += (dx[t0_i:t0_i+sample_window]-init_dx)**2 + (dy[t0_i:t0_i+sample_window]-init_dy)**2 + (dz[t0_i:t0_i+sample_window]-init_dz)**2
    
ds_sq_avg_sampled /= len(t0_iter_list)

# --------------------------------------------------------------------------------------------

fit_params, fit_cov = curve_fit(f, t_shortened, ds_sq_avg_sampled)
m, c = fit_params
err_m, err_c = np.sqrt(np.diag(fit_cov))

print('Fitted parameters:\nm = {:.4e}\nc = {:.4e}'.format(fit_params[0], fit_params[1]))
D = fit_params[0]/(6)
D_err = err_m/6

print('Diffusion coefficient:\nD = {:.4e} +/- {:.4e}'.format(D, D_err))
print('D = {:.4f} +/- {:.4f}'.format(D, D_err))

fitline = f(t_shortened, *fit_params)

# --------------------------------------------------------------------------------------------

# ---------------- TAU CALCULATION ----------------

tau = (D * (1E-9)**2) / (D_expected * (1E-6)**2)
print("Assuming\nD_experiment = {:.4f} micron^2/sec".format(D_expected))
print("We get")
print("tau = {:.4e} sec".format(tau))

tau = tau / (1E-6)

print("tau = {:.4f} microsecond".format(tau))

# --------------------------------------------------------------------------------------------

plt.figure(tight_layout=True, figsize=(5, 5))

plt.plot(t_shortened, ds_sq_avg_sampled, 'k-', label='Data')
plt.plot(t_shortened, fitline, 'r--', label='Fit')

plt.xlabel(r'$t/\tau$', fontsize=18)
plt.ylabel(r'$\langle \Delta s^2 \rangle$', fontsize=18)

plt.xlim(0, t_shortened[-1])
plt.ylim(bottom=0)

plt.title(r"$D_{{\mathrm{{CoM}}}} = {:.4f}$, $\tau = {:.4f}\,\mu s$".format(D, tau), fontsize=18)

plt.legend(fontsize=14)

plt.savefig('plots/D_com_sampled.pdf')

# --------------------------------------------------------------------------------------------