import numpy as np
import matplotlib.pyplot as plt
import rcparams

from scipy.signal import find_peaks
from scipy.optimize import brentq

s, corr_s = np.loadtxt('lp_data.dat', unpack=True)

zeros = []
for s_i in range(len(s) - 1):
    if corr_s[s_i] * corr_s[s_i + 1] < 0:
        zero = brentq(lambda x: corr_s[s_i] + (corr_s[s_i + 1] - corr_s[s_i]) * (x - s[s_i]) / (s[s_i + 1] - s[s_i]), s[s_i], s[s_i + 1])
        if zero <= max(s):
            zeros.append(zero)
            print("Zero found at s = {:.2f}".format(zero))
        
        
# Plot the correlation function and the zeros
plt.figure(figsize=(6, 6), constrained_layout=True)

plt.plot(s, corr_s, 'o', markersize=2, color='k', label='Correlation function', lw=0.5)
plt.axhline(0, color='k', lw=0.5, ls='--', alpha=0.5)
plt.axvline(0, color='k', lw=0.5, ls='--', alpha=0.5)
plt.scatter(zeros, [0] * len(zeros), color='r', label='Zeros', s=10)
plt.xlabel(r'$s$')
plt.ylabel("Correlation function")
plt.legend()
plt.savefig('zeros_distance.pdf', dpi=300)

# Use distances between zeros to calculate the radius of curvature

max_s = 150

R_estimates = []
for s_i in range(len(zeros)):
    for s_j in range(s_i, len(zeros)):
        if s_i != s_j and corr_s[s_i] < max_s and corr_s[s_j] < max_s:
            diff = np.abs(s_i - s_j)
            dist = np.abs(zeros[s_i] - zeros[s_j])
            R_estimate = dist / (diff * np.pi)
            R_estimates.append(R_estimate)
            
R_estimates = np.array(R_estimates)
R_mean = np.mean(R_estimates)
R_std = np.std(R_estimates)
print("Estimated radius of curvature: {:.2f} +/- {:.2f}".format(R_mean, R_std))