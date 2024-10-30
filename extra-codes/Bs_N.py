import numpy as np
import matplotlib.pyplot as plt
import rcparams

Ebind0 = -(4.8)
T = 1
beta = 1.0 / T
R0 = 100
Rb = 350
r_mono = 2.5
# Bs = 10000
N_list = np.arange(2, 100)
Bs_list = np.linspace(30000, 40000, 100)

normalize_lengths = 1

total_points_to_compute = len(Bs_list) * len(N_list)

Bs_N_matrix = np.zeros((len(Bs_list), len(N_list)))

points_computed = 0
for Bs_i, Bs in enumerate(Bs_list):
    # average_length_list = []
    for N_i, N in enumerate(N_list):
        n_list = np.arange(1, N+1)
        p_list = np.zeros_like(n_list, dtype=float)

        for n_i, n in enumerate(n_list):
            dEbend = r_mono * (n - 1) * Bs * ( (1/Rb) - (1/R0)  )**2
            dEbind = n * Ebind0
            p_list[n_i] = (N - n + 1) * np.exp(-beta * (dEbend + dEbind))

        Z = np.sum(p_list)

        p_list = p_list / Z

        average_n = np.sum(n_list * p_list)
        if normalize_lengths:
            average_n = average_n / N
        
        Bs_N_matrix[Bs_i, N_i] = average_n
        
        points_computed += 1
        percent_done = 100 * points_computed / total_points_to_compute
        if points_computed % 10000 == 0 or points_computed == total_points_to_compute:
            print("{} / {} | {:.2f}%".format(points_computed, total_points_to_compute, percent_done))

plt.figure()

if normalize_lengths:
    plt.imshow(Bs_N_matrix, aspect='auto', origin='lower', extent=[N_list[0], N_list[-1], Bs_list[0], Bs_list[-1]], cmap='magma', vmin=0, vmax=1)
else:
    plt.imshow(Bs_N_matrix, aspect='auto', origin='lower', extent=[N_list[0], N_list[-1], Bs_list[0], Bs_list[-1]], cmap='magma')

plt.colorbar()

plt.xlabel(r'$N$')
plt.ylabel(r'$B_s$')

plt.savefig('Bs_N_matrix.pdf')
    
    
    
    