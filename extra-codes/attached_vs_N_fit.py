import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import rcparams

n_input, n_attach_input = np.loadtxt('n_vs_attached_all.dat', unpack=True)
fraction_attached = n_attach_input / n_input

data_average = 0

counter = 0
for n, n_attach in zip(n_input, n_attach_input):
    if n > 10:
        data_average += n_attach
        counter += 1

data_average = data_average / counter

n_input = np.array(n_input)
n_attach_input = np.array(n_attach_input)
fraction_attached = np.array(fraction_attached)

# ------------------------------------------------------------------------

Rb = 350
R0 = 100
r_mono = 2.5
T = 1
kB = 1

def average_length(N_max, lp, Ebind0):
    n_list = np.arange(1, N_max+1)
    p_list = np.zeros_like(n_list, dtype=float)
    
    beta = 1.0 / (kB * T)

    for n_i, n in enumerate(n_list):
        dEbend = r_mono * (n - 1) * (lp * kB * T) * ( (1/Rb) - (1/R0)  )**2
        dEbind = n * Ebind0
        p_list[n_i] = (N_max - n + 1) * np.exp(-beta * (dEbend + dEbind))

    Z = np.sum(p_list)

    p_list = p_list / Z

    average_n = np.sum(n_list * p_list)

    return average_n

average_length_vec = np.vectorize(average_length)

# popt, pcov = curve_fit(lambda N_max, lp, Ebind0: average_length_vec(N_max, lp, Ebind0), n_input, n_attach_input, p0=[1200, -1.0], bounds=([0, -np.inf], [np.inf, 0]))

# lp, Ebind0 = popt
# err_lp, err_Ebind0 = np.sqrt(np.diag(pcov))
# print("lp = {:.2f} +/- {:.2f}".format(lp, err_lp))
# print("Ebind0 = {:.2f} +/- {:.2f}".format(Ebind0, err_Ebind0))

# average_length_fit = average_length_vec(n_input, lp, Ebind0)

n_attach_input = n_attach_input / n_input

plt.figure(figsize=(6, 4), dpi=300, constrained_layout=True)

plt.plot(n_input, n_attach_input, 'o', color='blue', label='Data')
# plt.plot(n_input, n_input, color='black', ls='--', lw=1)
plt.axhline(1, color='black', ls='--', lw=1, label=r'$\langle n \rangle/ N = 1$')
# plt.plot(n_input, average_length_fit, color='red', label=r'Fit:\\$l_p = {:.2f} \pm {:.2f}$\\$E_{{\mathrm{{bind}}}}^0 = {:.2f} \pm {:.2f}$'.format(lp, err_lp, Ebind0, err_Ebind0))

plt.xlabel(r'$N$', fontsize=16)
plt.ylabel(r'$\langle n \rangle / N$', fontsize=16)

plt.legend()

plt.savefig('attached_vs_N_fit.pdf')