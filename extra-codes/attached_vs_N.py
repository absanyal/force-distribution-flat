import numpy as np
import matplotlib.pyplot as plt
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

N_list = np.arange(4, 1000)
average_length_list = []

Ebind0 = -(0.1)
T = 1
kB = 1
beta = 1.0 / (kB * T)
R0 = 100
Rb = 350
r_mono = 2.5

lp = 1200

for N in N_list:

    n_list = np.arange(1, N+1)
    p_list = np.zeros_like(n_list, dtype=float)

    for n_i, n in enumerate(n_list):
        dEbend = r_mono * (n - 1) * (kB * T * lp) * ( (1/Rb) - (1/R0)  )**2
        dEbind = n * Ebind0
        p_list[n_i] = (N - n + 1) * np.exp(-beta * (dEbend + dEbind))

    Z = np.sum(p_list)

    p_list = p_list / Z

    average_n = np.sum(n_list * p_list)
    
    average_length_list.append(average_n)

    # plt.figure(figsize=(6, 4), dpi=300, constrained_layout=True)
    # plt.axvline(average_n, color='red', linestyle='--', label='Average: {:.2f}'.format(average_n))

    # plt.bar(n_list, p_list, width=1.0, color='blue')
    # plt.xlabel('Number of bound monomers')
    # plt.ylabel('Probability')

    # plt.legend()

    # plt.savefig('plot.pdf')

max_val = np.max(average_length_list)

plt.figure(figsize=(6, 4), dpi=300, constrained_layout=True)

plt.plot(N_list, average_length_list, color='black')


plt.xlabel('Number of monomers')
plt.ylabel('Average length attached')

# plt.xlim(N_list[0], N_list[-1])
# plt.xticks(np.arange(0, N_list[-1], 100))

# plt.ylim(bottom=0)

# plt.xscale('log')

# plt.axhline(max_val, color='red', linestyle='--', linewidth = 0.5, label='Max: {:.2f}'.format(max_val))

# plt.axhline(data_average, color='green', linestyle='--', linewidth = 0.5, label='Data max: {:.2f}'.format(data_average))

# plt.plot(n_input, n_attach_input, 'o', color='blue', label='Data')

# plt.plot(n_input, fraction_attached, 'o', color='blue', label='Data')

plt.title(r'$l_p = {:.2f}$ nm, $E_{{\mathrm{{bind}}}}^0 = {:.2f}$ $k_BT$'.format(lp, Ebind0))

# plt.legend()

plt.savefig('attached_length_vs_N.pdf')