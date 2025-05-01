import numpy as np
import matplotlib.pyplot as plt
import rcparams
from scipy.optimize import curve_fit

from numpy import arcsin, sqrt, pi

def delta_l(a, Rcell, R0, ks, ktheta):
    l_old = a * R0 / sqrt(4 * R0**2 - a**2)
    l_new = a * Rcell / sqrt(4 * Rcell**2 - a**2)
    return l_new - l_old

def percent_l_change(a, Rcell, R0, ks, ktheta):
    l_old = a * R0 / sqrt(4 * R0**2 - a**2)
    l_new = a * Rcell / sqrt(4 * Rcell**2 - a**2)
    return (l_new - l_old) / l_old

Rcell = 350
R0 = np.linspace(50, 1000, 1000)
a = 5

ks = 100
ktheta = 20

d_l = delta_l(a, Rcell, R0, ks, ktheta)
d_l_percent = percent_l_change(a, Rcell, R0, ks, ktheta) * 100

x_axis = np.zeros(len(R0))
for i in range(len(R0)):
    x_axis[i] = ( 1/R0[i] - 1/Rcell )**2
    
plt.figure(figsize=(6, 6))

plt.plot(x_axis, d_l_percent, 'o', label='Data', markersize=1)

plt.legend()

plt.show()