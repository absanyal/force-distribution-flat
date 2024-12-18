import numpy as np
import matplotlib.pyplot as plt

R_cell = 350

R0_list, Eb_list = np.loadtxt("measured_Ebend_vs_R0.dat", unpack=True)

x_list = ((1/R_cell) - (1/R0_list))**2

plt.figure()

plt.plot(x_list, Eb_list, 'o', label='Measured data', color='black')

plt.xlabel(r'$(1/R_{\rm cell} - 1/R_0)^2$')
plt.ylabel(r'$E_{\mathrm{bend}}\,(k_BT)$')

plt.legend()

plt.show()