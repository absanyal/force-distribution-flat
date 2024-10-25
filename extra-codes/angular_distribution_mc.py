import numpy as np
import matplotlib.pyplot as plt
import rcparams

from numpy import sin, cos, pi
from scipy.signal import find_peaks

Rc = 350
R0 = 100
r_mono = 2.5

T = 1
beta = 1.0 / T

n = 10

Bs = 1000
Lb = 2 * r_mono * n


N_iters = 1000000

theta = np.random.uniform(0, 2 * pi)

theta_list = []

fig = plt.figure(figsize=(10, 5), constrained_layout=True)
ax1 = fig.add_subplot(121, polar=True)
ax2 = fig.add_subplot(122)

for iter in range(1, N_iters+1):
    E_old = 0.5 * Bs * Lb * ((((cos(theta))**2) / Rc) - (1/R0))**2
    is_in_range = False
    while not is_in_range:
        # theta_new = theta + np.random.uniform(-pi/4, pi/4) * 0.01
        theta_new = np.random.uniform(0, 2 * pi)
        if 0 <= theta_new <= 2 * pi:
            is_in_range = True
    E_new = 0.5 * Bs * Lb * ((((cos(theta_new))**2) / Rc) - (1/R0))**2
    dE = E_new - E_old
    boltzmann = np.exp(-beta * dE)
    r = np.random.uniform(0, 1)
    if dE < 0 or r < boltzmann:
        theta = theta_new
    
    theta_list.append(theta)
    
    if iter % 100000 == 0:
        percent = (iter / N_iters) * 100
        print("Iter: {} / {} | {:.2f}%".format(iter, N_iters, percent))
    
# if iter % 5000 == 0:
ax1.clear()
ax1.hist(theta_list, bins='auto', density=True, color='b')

ax2.clear()
ax2.hist(theta_list, bins='auto', density=True, color='b')

plt.savefig("angular_distribution_mc.pdf")
