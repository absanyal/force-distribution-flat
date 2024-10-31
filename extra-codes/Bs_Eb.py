import numpy as np
import matplotlib.pyplot as plt
import rcparams

T = 1
beta = 1.0 / T
R0 = 100
Rb = 350
r_mono = 2.5
N = 2

tolerance = 0.01

normalize_lengths = 1

linker_multliplicity = 1


Bs_points = 100
Bs_min = 20000
Bs_max = 50000

Ebind0_points = 100
Ebind0_min = -10
Ebind0_max = 0

Bs_list = np.linspace(Bs_min, Bs_max, Bs_points)
Ebind0_list = np.linspace(Ebind0_min, Ebind0_max, Ebind0_points)


total_points_to_compute = len(Bs_list) * len(Ebind0_list)

Bs_Eb_matrix = np.zeros((len(Bs_list), len(Ebind0_list)), dtype=float)

points_computed = 0
for Bs_i, Bs in enumerate(Bs_list):
    # average_length_list = []
    for Ebind0_i, Ebind0 in enumerate(Ebind0_list):
        n_list = np.arange(1, N+1)
        p_list = np.zeros_like(n_list, dtype=float)

        for n_i, n in enumerate(n_list):
            dEbend = r_mono * (n - 1) * Bs * ( (1/Rb) - (1/R0)  )**2
            dEbind = n * Ebind0
            p_list[n_i] = (N - n + 1) * np.exp(-beta * (dEbend + dEbind))

        Z = np.sum(p_list)

        p_list = p_list / Z

        average_n = np.sum(n_list * p_list)
        average_n_unnormalized = average_n
        if normalize_lengths:
            average_n = average_n / N
            
        Bs_Eb_matrix[Bs_i, Ebind0_i] = average_n
        
        points_computed += 1
        percent_done = 100 * points_computed / total_points_to_compute
        if points_computed % 10000 == 0 or points_computed == total_points_to_compute:
            print("{} / {} | {:.2f}%".format(points_computed, total_points_to_compute, percent_done))

plt.figure()

if normalize_lengths:
    plt.imshow(Bs_Eb_matrix, aspect='auto', origin='lower', extent=[Ebind0_list[0], Ebind0_list[-1], Bs_list[0], Bs_list[-1]], cmap='magma')
else:
    plt.imshow(Bs_Eb_matrix, aspect='auto', origin='lower', extent=[Ebind0_list[0], Ebind0_list[-1], Bs_list[0], Bs_list[-1]], cmap='magma')

plt.colorbar()

plt.title(r'$N = {}$'.format(N))

plt.xlabel(r'$E_{\mathrm{bind}}^0$')
plt.ylabel(r'$B_s$')

plt.savefig('Bs_Eb_matrix.pdf')