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
        
n_input = n_input[1:]
n_attach_input = n_attach_input[1:]
fraction_attached = fraction_attached[1:]

n_input = np.array(n_input)
n_attach_input = np.array(n_attach_input)
fraction_attached = np.array(fraction_attached)

# ------------------------------------------------------------------------

Rb = 350
R0 = 100
r_mono = 2.5
T = 1

def average_length(N_max, Bs, Ebind0):
    n_list = np.arange(1, N_max+1)
    p_list = np.zeros_like(n_list, dtype=float)
    
    p_list = np.zeros_like(n_list, dtype=float)

    beta = 1.0 / T

    for n_i, n in enumerate(n_list):
        dEbend = r_mono * (n - 1) * Bs * ( (1/Rb) - (1/R0)  )**2
        dEbind = n * Ebind0
        p_list[n_i] = (N_max - n + 1) * np.exp(-beta * (dEbend + dEbind))

    Z = np.sum(p_list)

    p_list = p_list / Z

    average_n = np.sum(n_list * p_list)

    return average_n

average_length_vec = np.vectorize(average_length)

popt, pcov = curve_fit(average_length_vec, n_input, n_attach_input, p0=[5000, -1], maxfev=100000,    bounds=([0, -10], [100000, 0]))

average_length_fit = average_length_vec(n_input, *popt)

Bs, Ebind0 = popt
err_Bs, err_Ebind0 = np.sqrt(np.diag(pcov))
print("Bs = {:.2f} +/- {:.2f}".format(Bs, err_Bs))
print("Ebind0 = {:.2f} +/- {:.2f}".format(Ebind0, err_Ebind0))

plt.figure(figsize=(6, 4), dpi=300, constrained_layout=True)

plt.plot(n_input, n_attach_input, 'o', color='blue', label='Data')
plt.plot(n_input, average_length_fit, color='red', label='Fit: Bs = {:.2f}, Ebind0 = {:.2f}'.format(Bs, Ebind0))

plt.xlabel('Number of monomers')
plt.ylabel('Number of attached monomers')

plt.legend()

plt.savefig('attached_vs_N_fit.pdf')

test = average_length_vec(100, 38650, -4.8)
print(test)