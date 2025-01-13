import numpy as np
import matplotlib.pyplot as plt
import rcparams
from matplotlib.ticker import StrMethodFormatter

kBT0 = 310
sigma0 = 1

T = 310

epsilon = 800
sigma = 2.1
delta_E = 1.0

sigma = sigma / sigma0
T = T / kBT0
epsilon = epsilon / kBT0

e_req = -epsilon + delta_E


def LJ93(r, epsilon, sigma):
    return epsilon * ((2/15) * (sigma/r)**9 - (sigma/r)**3)

def LJ126(r, epsilon, sigma):
    return 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)

def LJ93_minus_e_req(r, epsilon, sigma, e_req):
    return LJ93(r, epsilon, sigma) - e_req

def LJ126_minus_e_req(r, epsilon, sigma, e_req):
    return LJ126(r, epsilon, sigma) - e_req

def p_LJ93(r, epsilon, sigma, T):
    beta = 1/T
    w_list = np.exp(-beta * LJ93(r, epsilon, sigma))
    Z = np.trapz(w_list, r)
    return w_list / Z

def p_LJ126(r, epsilon, sigma, T):
    beta = 1/T
    w_list = np.exp(-beta * LJ126(r, epsilon, sigma))
    Z = np.trapz(w_list, r)
    return w_list / Z

def find_roots(xlist, ylist):
    roots = []
    for i in range(len(ylist) - 1):
        if ylist[i] * ylist[i + 1] < 0:
            roots.append(xlist[i])
    return roots

def find_minima(xlist, ylist):
    minima = []
    for i in range(1, len(ylist) - 1):
        if ylist[i] < ylist[i - 1] and ylist[i] < ylist[i + 1]:
            minima.append(xlist[i])
    return minima
    
r = np.linspace(1, 4, 1000)
p_list = p_LJ93(r, epsilon, sigma, T)
V_list = LJ93(r, epsilon, sigma)

minimum = find_minima(r, V_list)

shifted_V_list = LJ93_minus_e_req(r, epsilon, sigma, e_req)

roots = find_roots(r, shifted_V_list)
root1, root2 = min(roots), max(roots)

p_in_allowed_area = np.trapz(p_list[(r >= root1) & (r <= root2)], r[(r >= root1) & (r <= root2)])

plt.figure(figsize=(6, 6))

plt.plot(r, p_list, label=r'$p(r)$', color='b', ls='--')
plt.plot(r, V_list, label=r'$V(r)$', color='r')

plt.axhline(y=0, color='k', linestyle='--', alpha=0.5, lw=0.5)
plt.axvline(x=0, color='k', linestyle='--', alpha=0.5, lw=0.5)

plt.axhline(y=e_req, color='r', linestyle='--', label=r'$E_{{\mathrm{{req}}}} = {:.4f} \,\mathrm{{kBT}}$'.format(e_req))

plt.axvline(x=minimum[0], color='r', linestyle='--', label=r'$r_{{0}} = {:.4f}\,\mathrm{{nm}}$'.format(minimum[0]), lw=0.5)

plt.axvline(x=root1, color='g', linestyle='--', label=r'$r_{{\mathrm{{min}}}} = {:.4f}\,\mathrm{{nm}}$'.format(root1))
plt.axvline(x=root2, color='g', linestyle='--', label=r'$r_{{\mathrm{{max}}}} = {:.4f}\,\mathrm{{nm}}$'.format(root2))

plt.axvspan(root1, root2, alpha=0.1, color='g', label=r'$p(r_{{\mathrm{{min}}}} \le r \le r_{{\mathrm{{max}}}})= {:.4f}$'.format(p_in_allowed_area))

plt.ylim(-epsilon - 1.0, max(p_list) + 1.0)
plt.xlim(min(r), max(r))

plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}')) # 2 decimal places

plt.xlabel(r'$r (\mathrm{nm})$')

plt.legend(fontsize=12)

plt.savefig('potential_and_probability.pdf')

######################################################################################

p_a = p_in_allowed_area

N = 20

# probability of n particles in the allowed area
p_n = np.zeros(N + 1)
p_n[0] = (1 - p_a)**N
for n in range(1, N + 1):
    p_n[n] = (1 - p_a)**(N - n) * p_a**n

normalization = np.sum(p_n)

p_n = p_n / normalization

average_n = np.sum(np.arange(N + 1) * p_n)

plt.figure(figsize=(5, 5))

plt.plot(range(N + 1), p_n, 'o-', color='k')

plt.axvline(x=average_n, color='r', linestyle='--', label=r'$\langle n \rangle = {:.2f}$'.format(average_n))

plt.xlabel(r'$n$')
plt.ylabel(r'$p(n)$')

plt.legend()

plt.savefig('probability_of_n_particles.pdf')