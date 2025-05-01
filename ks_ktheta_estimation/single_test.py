import numpy as np
import matplotlib.pyplot as plt
import rcparams
from scipy.optimize import curve_fit

from numpy import arcsin, sqrt, pi

def linear(x, m):
    return m * x

def dE_per_monomer(a, Rcell, R0, ks, ktheta):
    a1_old = a + ((a**2)) / sqrt(4 * R0**2 - a**2)
    a2_old = a - ((a**2)) / sqrt(4 * R0**2 - a**2)
    
    l_old = a * R0 / sqrt(4 * R0**2 - a**2)
    
    theta1_old = arcsin(a / (2 * l_old))
    theta2_old = pi - theta1_old
    
    a1_new = a + ((a**2)) / sqrt(4 * Rcell**2 - a**2)
    a2_new = a - ((a**2)) / sqrt(4 * Rcell**2 - a**2)
    
    l_new = a * Rcell / sqrt(4 * Rcell**2 - a**2)
    
    theta1_new = arcsin(a / (2 * l_new))
    theta2_new = pi - theta1_new
    
    d_a1 = a1_new - a1_old
    d_a2 = a2_new - a2_old
    d_theta1 = theta1_new - theta1_old
    d_theta2 = theta2_new - theta2_old
    
    dE = (ks * (d_a1**2 + d_a2**2) + ktheta * (d_theta1**2 + d_theta2**2))
    
    return dE

a = 5
Rcell = 350
R0 = np.arange(50, 501, 50)

ks = 100
ktheta = 20
    
dE = np.zeros(len(R0))

x_axis = np.zeros(len(R0))
for i in range(len(R0)):
    dE[i] = dE_per_monomer(a, Rcell, R0[i], ks, ktheta)
    x_axis[i] = ( 1/R0[i] - 1/Rcell )**2

popts, pcov = curve_fit(linear, x_axis, dE)
m = popts[0]
lp = m / (2 * (a/2))

print(("lp = {:.2f} nm").format(lp))

plt.figure(figsize=(6, 6))

plt.plot(x_axis, dE, 'o', label='Data')
plt.plot(x_axis, linear(x_axis, m), 'r-', label=r'Fit: $l_p = {:.2f}$ nm'.format(lp))

# plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

plt.xlabel(r'$(\frac{1}{R_0} - \frac{1}{R_{cell}})^2$ $({{\mathrm{{nm}}}}^{-2})$')
plt.ylabel(r'$\Delta E_{{\mathrm{{bend}}}}\,(k_B T)$')

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.savefig('single_test.pdf')