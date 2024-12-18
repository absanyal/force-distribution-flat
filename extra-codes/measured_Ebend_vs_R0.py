import numpy as np
import matplotlib.pyplot as plt
import rcparams
from scipy.optimize import curve_fit

def f(x, m, c):
    return m*x + c

R_cell = 350

r_mono = 2.5
n = 20

R0_list, Eb_list = np.loadtxt("measured_Ebend_vs_R0.dat", unpack=True)

x_list = ((1/R_cell) - (1/R0_list))**2

pop, pcov = curve_fit(f, x_list, Eb_list)

x_fit = np.linspace(min(x_list), max(x_list), 1000)
y_fit = f(x_fit, *pop)

m, c = pop

plt.figure()

plt.plot(x_list, Eb_list, 'o', label='Measured data', color='black')
# plt.plot(x_fit, y_fit, label=r'Fit: $E_{{\mathrm{{bend}}}} = ({:.2f})(1/R_{{\mathrm{{cell}}}} - 1/R_0)^2 + ({:.2f})$'.format(m, c), color='red', linestyle='--')
plt.plot(x_fit, y_fit, label="Fitting", color='red', linestyle='--')

plt.xlabel(r'$(1/R_{\rm cell} - 1/R_0)^2\,(\mathrm{nm}^{-2})$')
plt.ylabel(r'$E_{\mathrm{bend}}\,(k_BT)$')

# plt.xscale('log')
# plt.yscale('log')

plt.legend()

plt.savefig("measured_Ebend_vs_R0.pdf")

lp = m / (r_mono * (n-1))

print("Persistence length = {:.2f} nm".format(lp))

plt.clf()

plt.plot(R0_list, Eb_list, 'o', label='Measured data', color='black')

plt.xlabel(r'$R_0\,(\mathrm{nm})$')
plt.ylabel(r'$E_{\mathrm{bend}}\,(k_BT)$')

plt.legend()

plt.savefig("measured_Ebend_vs_R0_raw.pdf")

plt.clf()

x_list = ((1/R0_list))

plt.plot(x_list, Eb_list, 'o', label='Measured data', color='black')

plt.xlabel(r'$1/R_0\,(\mathrm{nm}^{-1})$')
plt.ylabel(r'$E_{\mathrm{bend}}\,(k_BT)$')

plt.legend()

plt.savefig("measured_Ebend_vs_R0_raw_inverse.pdf")