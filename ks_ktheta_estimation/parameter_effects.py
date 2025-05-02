import numpy as np
import matplotlib.pyplot as plt
import rcparams
from scipy.optimize import curve_fit

from numpy import arcsin, sqrt, pi

colorlist = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta']
markerstyle_list = ['o', 's', '^', 'D', 'x', 'v', 'p', '*', '<', '>', 'h']

def linear(x, m):
    return m * x

def dE_per_monomer(a, Rcell, R0, ks, ktheta, kl):
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
    
    d_l = 2*(l_new - l_old)
    
    dE = (ks * (d_a1**2 + d_a2**2) + ktheta * (d_theta1**2 + d_theta2**2)) + 2 * kl * (d_l**2)
    
    return dE

def dE_full_filament(N, a, Rcell, R0, ks, ktheta, kl):
    num_a1 = 2 * N
    num_a2 = 2 * N
    num_l = 2 * (N+1)
    num_theta1 = 2 * N
    num_theta2 = 2 * N
    
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
    
    d_a1 = (a1_new - a1_old) * num_a1
    d_a2 = (a2_new - a2_old) * num_a2
    d_theta1 = (theta1_new - theta1_old) * num_theta1
    d_theta2 = (theta2_new - theta2_old) * num_theta2
    
    d_l = (2*(l_new - l_old)) * num_l
    
    dE_a1 = 0.5 * ks * (d_a1**2) * num_a1
    dE_a2 = 0.5 * ks * (d_a2**2) * num_a2
    dE_theta1 = 0.5 * ktheta * (d_theta1**2) * num_theta1
    dE_theta2 = 0.5 * ktheta * (d_theta2**2) * num_theta2
    dE_l = 0.5 * kl * (d_l**2) * num_l
    
    dE = dE_a1 + dE_a2 + dE_theta1 + dE_theta2 + dE_l
    
    return dE
    
    

Rcell = 350
# R0 = np.linspace(50, 1000, 1000)
R0 = np.array([100, 150, 200, 250, 300])
a = 5
N = 20

ks_list = np.array([20, 50, 100])
# ktheta_list = np.arange(1, 201, 50)
ktheta_list = np.array([10, 20])
kl = 100

x_axis = np.zeros((len(R0)))
for i in range(len(R0)):
            x_axis[i] = ( 1/R0[i] - 1/Rcell )**2


points_to_calculate = len(ks_list) * len(ktheta_list)
print("Total points to calculate: {}".format(points_to_calculate))
step_size = int(0.1 * points_to_calculate)

plt.figure(figsize=(6, 6))

# points_calculated = 0
for ks_i, ks in enumerate(ks_list):
    for ktheta_i, ktheta in enumerate(ktheta_list):
        
        dE = np.zeros(len(R0))
        
        for i in range(len(R0)):
            # dE[i] = dE_full_filament(N, a, Rcell, R0[i], ks, ktheta, kl)
            dE[i] = dE_per_monomer(a, Rcell, R0[i], ks, ktheta, kl)
        
        popts, pcov = curve_fit(linear, x_axis, dE)
        m = popts[0]
        
        # If using dE_per_monomer
        lp = m / ((a/2))
        
        # If using dE_full_filament
        # lp = m / ((a/2) * (N))
        print("ks = {}, ktheta = {} | lp = {:.2f} nm".format(ks, ktheta, lp))
        
        thiscolor = colorlist[ks_i % len(colorlist)]
        thismarker = "{}-".format(markerstyle_list[ktheta_i % len(markerstyle_list)])
        
        plt.plot(x_axis, dE, thismarker, label=r'$k_s = {}$, $k_\theta = {}$, $l_p = {:.2f}$'.format(ks, ktheta, lp),
                 markersize=3, lw=0.5, color=thiscolor)
        
        # points_calculated += 1
        # percentage = (points_calculated / points_to_calculate) * 100
        # if points_calculated % step_size == 0 or points_calculated == points_to_calculate and points_to_calculate >= 10:
        #     print("{} / {} | {:.2f} %".format(points_calculated, points_to_calculate, percentage))
            
plt.xlabel(r'$(\frac{1}{R_0} - \frac{1}{R_{cell}})^2$')
plt.ylabel(r'$\Delta E$ (kT)')

# plt.xscale('log')
# plt.yscale('log')

plt.legend(fontsize=8)

plt.savefig('ks_ktheta_effects.pdf')