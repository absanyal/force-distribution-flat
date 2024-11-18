import numpy as np
import matplotlib.pyplot as plt
import rcparams

ktheta, lp, lp_err, e2e, e2e_err = np.loadtxt('ktheta_vs_lp.dat', unpack=True)

fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

ax.errorbar(ktheta, lp, yerr=lp_err, fmt='none', color='r', elinewidth=1, capsize=2)

ax.plot(ktheta, lp, color='k', label='Data', marker='o', markersize=4)

ax.set_xlabel(r'$k_{\theta}\,(k_BT)$')
ax.set_ylabel(r'$l_p$ (nm)')

plt.savefig('ktheta_vs_lp.pdf')

fig, ax = plt.subplots(constrained_layout=True, figsize=(6, 4))

ax.errorbar(ktheta, e2e, yerr=e2e_err, fmt='none', color='r', elinewidth=1, capsize=2)

ax.plot(ktheta, e2e, color='k', label='Data', marker='o', markersize=4)

ax.set_xlabel(r'$k_{\theta}\,(k_BT)$')
ax.set_ylabel(r'$\langle R_{{e-e}} \rangle$ (nm)')

plt.savefig('ktheta_vs_e2e.pdf')