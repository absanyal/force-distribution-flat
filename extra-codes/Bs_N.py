import numpy as np
import matplotlib.pyplot as plt
import rcparams

Ebind0 = -1.0
kB = 1
T = 1
beta = 1.0 / (kB * T)
R0 = 100
Rb = 350
r_mono = 2.5

N_min = 2
N_max = 30

N_list = np.arange(N_min, N_max+1)

lp_points = 100
lp_min = 1000
lp_max = 15000

lp_list = np.linspace(lp_min, lp_max, lp_points)

normalize_lengths = 1

total_points_to_compute = len(lp_list) * len(N_list)

lp_N_matrix = np.zeros((len(lp_list), len(N_list)))

points_computed = 0
for lp_i, lp in enumerate(lp_list):
    # average_length_list = []
    for N_i, N in enumerate(N_list):
        n_list = np.arange(1, N+1)
        p_list = np.zeros_like(n_list, dtype=float)

        for n_i, n in enumerate(n_list):
            dEbend = r_mono * (n - 1) * (kB * T *lp) * ( (1/Rb) - (1/R0)  )**2
            dEbind = n * Ebind0
            p_list[n_i] = (N - n + 1) * np.exp(-beta * (dEbend + dEbind))

        Z = np.sum(p_list)

        p_list = p_list / Z

        average_n = np.sum(n_list * p_list)
        if normalize_lengths:
            average_n = average_n / N
        
        lp_N_matrix[lp_i, N_i] = average_n
        
        points_computed += 1
        percent_done = 100 * points_computed / total_points_to_compute
        if points_computed % 10000 == 0 or points_computed == total_points_to_compute:
            print("{} / {} | {:.2f}%".format(points_computed, total_points_to_compute, percent_done))

plt.figure()

if normalize_lengths:
    plt.imshow(lp_N_matrix, aspect='auto', origin='lower', extent=[N_list[0], N_list[-1], lp_list[0], lp_list[-1]], cmap='magma', vmin=0, vmax=1)
else:
    plt.imshow(lp_N_matrix, aspect='auto', origin='lower', extent=[N_list[0], N_list[-1], lp_list[0], lp_list[-1]], cmap='magma')

plt.colorbar()

plt.title(r'$E_{{\mathrm{{bind}}}}^0 = {:.2f}$'.format(Ebind0))

plt.xlabel(r'$N$')
plt.ylabel(r'$l_p$ (nm)')

plt.savefig('Bs_N_matrix.pdf')
    
    
    
    