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

ks_list = np.linspace(10, 200, 100)
ktheta_list = np.linspace(10, 200, 100)

lp_matrix = np.zeros((len(ks_list), len(ktheta_list)))

total_values = len(ks_list) * len(ktheta_list)
step = int(0.1 * total_values)

evaluated_values = 0
for ktheta_i, ktheta in enumerate(ktheta_list):
    for ks_i, ks in enumerate(ks_list):
    
        dE = np.zeros(len(R0))

        x_axis = np.zeros(len(R0))
        for i in range(len(R0)):
            dE[i] = dE_per_monomer(a, Rcell, R0[i], ks, ktheta)
            x_axis[i] = ( 1/R0[i] - 1/Rcell )**2

        popts, pcov = curve_fit(linear, x_axis, dE)
        m = popts[0]
        lp = m / (2 * (a/2))
        
        lp_matrix[ks_i, ktheta_i] = lp
        
        evaluated_values += 1
        if evaluated_values % step == 0 or evaluated_values == total_values:
            percentage = evaluated_values / total_values * 100
            print("{} / {} | {:.2f}%".format(evaluated_values, total_values, percentage))
        

plt.figure(figsize=(8, 6))

plt.imshow(lp_matrix, aspect='auto', extent=[ktheta_list[0], ktheta_list[-1], ks_list[0], ks_list[-1]], origin='lower')

plt.colorbar(label=r'$l_p$ (nm)')

plt.xlabel(r'$k_{\theta}\,(k_B T)$')
plt.ylabel(r'$k_{s}\,(k_B T)$')

plt.savefig('persistence_length.pdf')