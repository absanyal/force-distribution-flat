from sys import argv
import numpy as np
from analysis_modules.cylindermath import cylinder, moving_average, distance_from_surface
from analysis_modules.create_cell import make_box
import matplotlib.pyplot as plt
import modules.rcparams

r0 = 1.8028
r_min = 1.6326
r_max = 2.1251

r_list_file, p_list_file = np.loadtxt('extra-codes/expected_pr.dat', unpack=True)

s_list = np.loadtxt('r_histogram.dat')

average_r = np.mean(s_list)

plt.figure(figsize=(5, 5), constrained_layout=True)

plt.hist(s_list, bins='auto', density=True, label=r'$\Delta\,s_{\mathrm{linker}}$ from LAMMPS', color='b', edgecolor='b', histtype='step', rwidth=1, lw=2.0)

plt.axvline(x=r0, color='r', linestyle='--', label=r'$r_0 = {:.4f}$ nm'.format(r0), lw=1.0)
plt.axvline(x=r_min, color='g', linestyle='--', label=r'$r_{{\mathrm{{min}}}} = {:.4f}$ nm'.format(r_min), lw=1.0)
plt.axvline(x=r_max, color='g', linestyle='--', label=r'$r_{{\mathrm{{max}}}} = {:.4f}$ nm'.format(r_max), lw=1.0)
plt.axvline(x=average_r, color='b', linestyle='--', label=r'$\langle\Delta\,s_{{\mathrm{{linker}}}}\rangle = {:.4f}$ nm'.format(average_r), lw=1.0)

plt.plot(r_list_file, p_list_file, '-', color='k', label='Expected distribution')

plt.xlabel(r'$\Delta\,s_{{\mathrm{{linker}}}}$ (nm)')
plt.ylabel('Probability density')

plt.xlim(left = 1, right = max(r_list_file))

plt.legend(fontsize=10)

plt.savefig('plots/r_hist.pdf')