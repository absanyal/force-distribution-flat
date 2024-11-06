import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import rcparams

def f(x, m, c):
    return m*x + c

gap, p = np.loadtxt('gap_vs_p_attach.dat', unpack=True)

gap = gap - 1

popt, pcov = curve_fit(f, gap, p)

m, c = popt
err_m, err_c = np.sqrt(np.diag(pcov))

print(r'm = {:.2f} +/- {:.2f}'.format(m, err_m))
print(r'c = {:.2f} +/- {:.2f}'.format(c, err_c))

test_gap = np.linspace(gap.min(), gap.max(), 100)
test_p = f(test_gap, m, c)

fig, ax = plt.subplots(constrained_layout=True)

ax.plot(gap, p, 'o-', color='k')

ax.plot(test_gap, test_p, '--', color='r', label='Fit: $y = {:.2f}x + {:.2f}$'.format(m, c))

ax.set_xlabel('Gap between linkers')
ax.set_ylabel('Probability of attachment')

ax.legend()

plt.savefig('gap_vs_p_attach.pdf')