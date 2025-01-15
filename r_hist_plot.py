from sys import argv
import numpy as np
from analysis_modules.cylindermath import cylinder, moving_average, distance_from_surface
from analysis_modules.create_cell import make_box
import matplotlib.pyplot as plt
import modules.rcparams

r0 = 1.8028
r_min = 1.6326
r_max = 2.1251

s_list = np.loadtxt('r_histogram.dat')

average_r = np.mean(s_list)

plt.figure(figsize=(5, 5))

plt.hist(s_list, bins=100, density=True, label=r'$\Delta\,s_{\mathrm{linker}}$', color='b', alpha=0.5, edgecolor='None', histtype='stepfilled', rwidth=0.8)

plt.axvline(x=r0, color='r', linestyle='--', label=r'$r_0 = {:.4f}$ nm'.format(r0))
plt.axvline(x=r_min, color='g', linestyle='--', label=r'$r_{{\mathrm{{min}}}} = {:.4f}$ nm'.format(r_min))
plt.axvline(x=r_max, color='g', linestyle='--', label=r'$r_{{\mathrm{{max}}}} = {:.4f}$ nm'.format(r_max))
plt.axvline(x=average_r, color='b', linestyle='--', label=r'$\langle\Delta\,s_{{\mathrm{{linker}}}}\rangle = {:.4f}$ nm'.format(average_r))

plt.xlabel(r'$\Delta\,s_{{\mathrm{{linker}}}}$ (nm)')
plt.ylabel('Probability density')

plt.legend()

plt.savefig('plots/r_hist.pdf')