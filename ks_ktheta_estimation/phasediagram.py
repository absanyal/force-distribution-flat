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

Rcell = 350
R0 = np.linspace(50, 1000, 1000)
a = 5

ks_list = np.arange(1, 101, 1)
ktheta_list = np.arange(1, 101, 1)

x_axis = np.zeros((len(R0)))
for i in range(len(R0)):
            x_axis[i] = ( 1/R0[i] - 1/Rcell )**2

ks_ktheta_phase_diagram = np.zeros((len(ks_list), len(ktheta_list)))


points_to_calculate = len(ks_list) * len(ktheta_list)
print("Total points to calculate: {}".format(points_to_calculate))
step_size = int(0.1 * points_to_calculate)

points_calculated = 0
for ks_i, ks in enumerate(ks_list):
    for ktheta_i, ktheta in enumerate(ktheta_list):
        
        dE = np.zeros(len(R0))
        
        for i in range(len(R0)):
            dE[i] = dE_per_monomer(a, Rcell, R0[i], ks, ktheta)
        
        popts, pcov = curve_fit(linear, x_axis, dE)
        m = popts[0]
        lp = m / (2 * (a/2))
        # print("ks = {}, ktheta = {} | lp = {:.2f} nm".format(ks, ktheta, lp))
        
        ks_ktheta_phase_diagram[ks_i, ktheta_i] = lp
        
        points_calculated += 1
        percentage = (points_calculated / points_to_calculate) * 100
        
        if points_calculated % step_size == 0 or points_calculated == points_to_calculate:
            print("{} / {} | {:.2f} %".format(points_calculated, points_to_calculate, percentage))
        

plt.figure(figsize=(8, 6))

plt.imshow(ks_ktheta_phase_diagram, aspect='auto', cmap='viridis', origin='lower', extent=[ktheta_list[0], ktheta_list[-1] + 1, ks_list[0], ks_list[-1] + 1])

# plt.imshow(ks_ktheta_phase_diagram, aspect='auto', cmap='viridis', origin='lower')

plt.xlabel(r'$k_\theta$ $(k_B T)$')
plt.ylabel(r'$k_s$ $(k_B T)$')

plt.colorbar(label=r'$l_p$ (nm)', orientation='vertical', pad=0.02)

plt.savefig('persistence_length_phase_diagram.pdf')

with open('ks_ktheta_phase_diagram.txt', 'w') as f:
    f.write("# ks ktheta lp\n")
    for ks_i, ks in enumerate(ks_list):
        for ktheta_i, ktheta in enumerate(ktheta_list):
            f.write("{} {} {}\n".format(ks, ktheta, ks_ktheta_phase_diagram[ks_i, ktheta_i]))