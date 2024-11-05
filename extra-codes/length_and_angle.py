import numpy as np
import matplotlib.pyplot as plt
import rcparams
from numpy import sin, cos, pi

N_list = np.arange(1, 100)
# theta_list = [0, pi/6, pi/3, pi/2, 2 * pi/3, 5 * pi/6, pi]
# theta_label_list = ['0', r'$\frac{\pi}{6}$', r'$\frac{\pi}{3}$', r'$\frac{\pi}{2}$', r'$\frac{2\pi}{3}$', r'$\frac{5\pi}{6}$', r'$\pi$']

theta_list = []
theta_label_list = []

theta_increment = 30
theta_min = 0
theta_max = 90

linker_multiplicity = 1

Ebind0_full = -4.8
Ebind0 = Ebind0_full / linker_multiplicity
T = 1
beta = 1.0 / T
R0 = 100
Rb = 350
r_mono = 2.5
Bs = 38500

for theta_i in range(theta_min, theta_max+1, theta_increment):
    theta_label_list.append(r'${}^\circ$'.format(theta_i))
    
    theta = np.deg2rad(theta_i)
    theta_list.append(theta)

plt.figure(figsize=(6, 4), dpi=300, constrained_layout=True)

for theta_i, theta in enumerate(theta_list):
    average_length_list = []
    for N_i, N in enumerate(N_list):

        n_list = np.arange(1, N+1)
        p_list = np.zeros_like(n_list, dtype=float)

        for n_i, n in enumerate(n_list):
            dEbend = r_mono * (n - 1) * Bs * ( ( ((cos(theta)**2)) /Rb) - (1/R0)  )**2
            dEbind = n * Ebind0
            p_list[n_i] = (N - n + 1) * np.exp(-beta * (dEbend + dEbind))

        Z = np.sum(p_list)

        p_list = p_list / Z

        average_n = np.sum(n_list * p_list) 
        
        average_length_list.append(average_n)

    max_val = np.max(average_length_list)

    # plt.figure(figsize=(6, 4), dpi=300, constrained_layout=True)
    plt.plot(N_list, average_length_list, label=r'{}, $n_{{\mathrm{{max}}}} = {:.2f}$'.format(theta_label_list[theta_i], max_val), linewidth=1.0)
    plt.xlabel('Number of monomers')
    plt.ylabel('Average length attached')

    # plt.xlim(N_list[0], N_list[-1])
    # plt.xticks(np.arange(0, N_list[-1], 100))
    
    plt.xscale('log')
    # plt.yscale('log')

    # plt.ylim(bottom=0)

    # plt.axhline(max_val, color='red', linestyle='--', linewidth = 0.5, label='Max: {:.2f}'.format(max_val))

# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
legend = plt.legend(loc='best')
legend.get_frame().set_alpha(0.2)

plt.title(r'$B_s = {:.2f}, E_{{\mathrm{{bind}}}}^0 = {:.2f}$'.format(Bs, Ebind0_full))

plt.savefig('length_and_angle.pdf')