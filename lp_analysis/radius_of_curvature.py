import numpy as np
import matplotlib.pyplot as plt
import rcparams

from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def correlation(s, lp, R):
    """
    Calculate the correlation function.
    """
    return np.exp(-s / lp) * np.cos(s / R)

s, corr_s = np.loadtxt('lp_data.dat', unpack=True)

max_s_prev, lp_prev, lp_log_prev = np.loadtxt('lp_envelope.dat', unpack=True)

L = len(s) + 1
print("L: {}".format(L))

# max_s = max(s)
max_s = 100

print("Max s: {}".format(max_s))

s = s[s < max_s]
max_index = len(s)
corr_s = corr_s[:max_index]



# lp_guess = 96.44
lp_guess = lp_log_prev
R_guess = 25

pop, pcov = curve_fit(correlation, s, corr_s, p0=[lp_guess, R_guess])


lp = pop[0]
lp_err = np.sqrt(np.diag(pcov))[0]
R = pop[1]
R_err = np.sqrt(np.diag(pcov))[1]

print("lp: {:.2f} +/- {:.2f}".format(lp, lp_err))
print("R: {:.2f} +/- {:.2f}".format(R, R_err))

y_fit = correlation(s, lp, R)

plt.figure(figsize=(6, 6), constrained_layout=True)

plt.plot(s, corr_s, 'o', markersize=2, color='k', label='Correlation function', lw=0.5)

plt.plot(s, y_fit, 'r-', label=r'Fit: $lp = {:.2f} \pm {:.2f}$, $R = {:.2f} \pm {:.2f}$'.format(lp, lp_err, R, R_err), lw=0.5)

plt.xlabel(r'$s$')
plt.ylabel("Correlation function")


plt.axhline(0, color='k', lw=0.5, ls='--', alpha=0.5)
plt.axvline(0, color='k', lw=0.5, ls='--', alpha=0.5)

plt.legend(loc='best', fontsize=8)

plt.savefig('lp_R_correlation.pdf', dpi=300)