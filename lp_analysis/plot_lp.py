import numpy as np
import matplotlib.pyplot as plt
import rcparams

from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def envelope(s, lp):
    """
    Calculate the envelope of the correlation function.
    """
    return np.exp(-s / lp)

def linear_envelope(s, lp):
    """
    Linear fit function for the log of the envelope
    """
    return -s / lp

s, corr_s = np.loadtxt('lp_data.dat', unpack=True)

L = len(s) + 1
print("L: {}".format(L))

# max_s = 60
# max_s = int(s[-1] * 0.5)
max_s = max(s)

print("Max s: {}".format(max_s))

s = s[s < max_s]
max_index = len(s)
corr_s = corr_s[:max_index]

peaks = find_peaks(corr_s)
peak_indices = peaks[0]

peak_indices = np.insert(peak_indices, 0, 0)

if len(peak_indices) > 1:
    print("Fitting using {} peaks".format(len(peak_indices)))

    peak_heights = [corr_s[i] for i in peak_indices]
    peak_heights = np.array(peak_heights)

    log_peaks = np.log(peak_heights)
else:
    print("Peaks found: {}".format(len(peak_indices)))
    print("Not enough peaks to fit the envelope.")

# Fit the envelope to the correlation function
if len(peak_indices) > 1:
    popt, pcov = curve_fit(envelope, s[peak_indices], peak_heights)

    lp = popt[0]
    lp_err = np.sqrt(np.diag(pcov))[0]
    print("lp: {:.2f} +/- {:.2f}".format(lp, lp_err))

    y_fit = envelope(s, lp)

    # Fit the linear envelope to the log of the correlation function
    popt_log, pcov_log = curve_fit(linear_envelope, s[peak_indices], log_peaks)
    lp_log = popt_log[0]
    lp_log_err = np.sqrt(np.diag(pcov_log))[0]
    print("lp_log: {:.2f} +/- {:.2f}".format(lp_log, lp_log_err))

    y_fit_log = np.exp(linear_envelope(s, lp_log))

#############################################################################################

if len(peak_indices) > 1:
    L_by_lp = L / lp
    print("L/lp: {:.2f}".format(L_by_lp))

    L_by_lp_log = L / lp_log
    print("L/lp_log: {:.2f}".format(L_by_lp_log))

###########################################################################################

plt.figure(figsize=(6, 6), constrained_layout=True)

plt.plot(s, corr_s, 'o-', markersize=2, color='k', label='Correlation function', lw=0.5)

if len(peak_indices) > 1:
    for i in peak_indices:
        plt.scatter(s[i], corr_s[i], color='red', s=50, marker='x')

    plt.plot(s, y_fit, 'b--', label=r'$l_p$ (exponential fit) = ${:.2f}$ $\pm$ ${:.2f}$\\ $L/l_p = {:.2f}$'.format(lp, lp_err, L_by_lp))

    plt.plot(s, y_fit_log, 'r--', label=r'$l_p$ (linear fit) = ${:.2f}$ $\pm$ ${:.2f}$\\ $L/l_p = {:.2f}$'.format(lp_log, lp_log_err, L_by_lp_log))

plt.xlabel(r'$s$')
# plt.ylabel(r'$C(s)$')
# plt.ylabel(r'$ \frac{1}{N-1} \sum \limits_{i} \langle \vec{\tau}_{i} \cdot \vec{\tau}_{i+s} \rangle$')
plt.ylabel('Correlation function')

plt.axhline(0, color='k', lw=0.5, ls='--', alpha=0.5)
plt.axvline(0, color='k', lw=0.5, ls='--', alpha=0.5)

# plt.yscale('log')
# plt.ylim(bottom=1e-4)

plt.legend()

plt.savefig('persistence_linkers.pdf', bbox_inches='tight')

with open('lp_envelope.dat', 'w') as f:
    f.write("# max_s lp lp_log\n")
    f.write("{} {} {}\n".format(max_s, lp, lp_log))