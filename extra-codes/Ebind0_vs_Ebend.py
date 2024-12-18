import numpy as np
import matplotlib.pyplot as plt
import rcparams
from scipy.optimize import curve_fit

def f(x, m, c):
    return m*x + c

R_cell = 350

r_mono = 2.5
n = 20

kBT0 = 310

Ebind0, Ebend = np.loadtxt("Ebind0_vs_Ebend.dat", unpack=True)

Ebind0 = Ebind0/kBT0

min_Ebend = np.min(Ebend)
max_Ebend = np.max(Ebend)
avg_Ebend = np.mean(Ebend)

plt.figure()

plt.plot(Ebind0, Ebend, 'o', label='Measured data', color='black')

plt.xlabel(r'$E_{\mathrm{bind}}^0\,(k_BT)$')
plt.ylabel(r'$E_{\mathrm{bend}}\,(k_BT)$')

plt.axhline(avg_Ebend, color='green', label=r'$\langle E_{{\mathrm{{bend}}}}\rangle = {:.4f}\,k_BT$'.format(avg_Ebend), linestyle='--')
plt.axhline(min_Ebend, color='red', label=r'$\min(E_{{\mathrm{{bend}}}}) = {:.4f}\,k_BT$'.format(min_Ebend), linestyle='--')
plt.axhline(max_Ebend, color='blue', label=r'$\max(E_{{\mathrm{{bend}}}}) = {:.4f}\,k_BT$'.format(max_Ebend), linestyle='--')

plt.legend()

plt.savefig("Ebind0_vs_Ebend.pdf")