import numpy as np
import matplotlib.pyplot as plt
import rcparams
from scipy.optimize import curve_fit
from numpy import arcsin, sqrt, pi

colorlist = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta', 'olive', 'lime', 'teal', 'lavender', 'salmon', 'turquoise', 'tan', 'gold', 'orchid', 'skyblue', 'lightgreen', 'lightblue', 'lightcoral', 'lightpink', 'lightgray', 'lightcyan', 'lightyellow', 'lightpurple', 'lightbrown', 'darkred', 'darkblue', 'darkgreen', 'darkorange', 'darkpurple', 'darkbrown', 'darkpink', 'darkgray', 'darkcyan', 'darkmagenta', 'darkolive', 'darklime', 'darkteal', 'darklavender', 'darksalmon', 'darkturquoise', 'darktan', 'darkgold', 'darkorchid', 'darkskyblue', 'darklightgreen', 'darklightblue', 'darklightcoral', 'darklightpink', 'darklightgray', 'darklightcyan', 'darklightyellow', 'darklightpurple', 'darklightbrown']
markerstyle_list = ['o', 's', '^', 'D', 'x', 'v', 'p', '*', '<', '>', 'h']


def f(x, m):
    return m*x


def read_data(filename):
    with open(filename, 'r') as f:
        datasets = []
        groupname_list = []
        current_dataset_x = []
        current_dataset_y = []
        for line in f:
            line = line.strip()
            if len(line) == 0:
                if len(current_dataset_x) > 0 and len(current_dataset_y) > 0:
                    datasets.append([current_dataset_x, current_dataset_y])
                    current_dataset_x = []
                    current_dataset_y = []
            elif line[0] == '!':
                groupname = line[1:]
                groupnames = groupname.split()
                groupnames = np.array(groupnames, dtype=float)
                groupname_list.append(groupnames)
            elif line[0] == '#':
                continue
            else:
                xval, yval = line.split()
                xval = float(xval)
                Ebend = float(yval)
                current_dataset_x.append(xval)
                current_dataset_y.append(Ebend)
        datasets.append([current_dataset_x, current_dataset_y])

        groupname_list = np.array(groupname_list)
        groupname_list = groupname_list.transpose()

    return groupname_list, datasets

def linear(x, m):
    return m * x

def dE_per_monomer(a, Rcell, R0, ks, ktheta, kl):
    a1_old = a + ((a**2)) / sqrt(4 * R0**2 - a**2)
    a2_old = a - ((a**2)) / sqrt(4 * R0**2 - a**2)
    
    l_old = a * R0 / sqrt(4 * R0**2 - a**2)
    
    theta1_old = arcsin(a / (2 * l_old))
    theta2_old = pi - theta1_old
    
    a1_new = a + ((a**2)) / sqrt(4 * Rcell**2 - a**2)
    a2_new = a - ((a**2)) / sqrt(4 * Rcell**2 - a**2)
    
    l_new = a * Rcell / sqrt(4 * Rcell**2 - a**2)
    
    theta1_new = arcsin(a / (2 * l_new))
    theta2_new = pi - theta1_new
    
    d_a1 = a1_new - a1_old
    d_a2 = a2_new - a2_old
    d_theta1 = theta1_new - theta1_old
    d_theta2 = theta2_new - theta2_old
    
    d_l = 2*(l_new - l_old)
    
    dE = (ks * (d_a1**2 + d_a2**2) + ktheta * (d_theta1**2 + d_theta2**2)) + 2 * kl * (d_l**2)
    
    return dE

headers, datasets = read_data("measured_Ebend_vs_R0.dat")

ks_list, ktheta_list = headers

ks_banlist = []
ks_banlist = np.array(ks_banlist)

R_cell = 350
kBT0 = 1

kl = 100
a = 5

r_mono = a / 2
n = 20

R0_list = np.array([100, 150, 200, 250, 300])
R0_list = np.array(R0_list)

x_axis = np.zeros((len(R0_list)))
for i in range(len(R0_list)):
            x_axis[i] = ( 1/R0_list[i] - 1/R_cell )**2

plt.figure(figsize=(10, 5), constrained_layout=True)

for ks_i in range(len(ks_list)):

    ks = ks_list[ks_i]
    ktheta = ktheta_list[ks_i]
    thiscolor = colorlist[ks_i % len(colorlist)]
    thismarker = "{}".format(markerstyle_list[ks_i % len(markerstyle_list)])

    if ks not in ks_banlist:

        ks = ks / kBT0
        ktheta = ktheta / kBT0

        print("ks = {:.2f} kBT, ktheta = {:.2f} kBT".format(ks, ktheta))

        R0_list, Eb_list = datasets[ks_i]

        R0_list = np.array(R0_list)
        Eb_list = np.array(Eb_list)

        x_list = ((1/R_cell) - (1/R0_list))**2

        pop, pcov = curve_fit(f, x_list, Eb_list)

        x_fit = np.linspace(min(x_list), max(x_list), 1000)
        y_fit = f(x_fit, *pop)

        m = pop[0]
        err_m = np.sqrt(np.diag(pcov))[0]
        
        lp = m / (r_mono * (n-1))
        lp_err = err_m / (r_mono * (n-1))

        plt.plot(x_list, Eb_list, thismarker, color=thiscolor, markersize = 3)
        plt.plot(x_fit, y_fit, color=thiscolor, linestyle='--', lw = 1, label=r'$k_s = {:.2f}\,k_BT,\, k_{{\theta}} = {:.2f}\,k_BT$, $l_p = {:.2f} \pm {:.2f}\,\mathrm{{nm}}$'.format(ks, ktheta, lp, lp_err))

        # print("Measured slope = {:.2f}".format(m))
        print("lp (fitting)= {:.2f} +/- {:.2f} nm".format(lp, lp_err))
        
        ################################################################################
        
        # Calculate the persistence length using the theoretical formula
        
        dE_list = np.zeros(len(R0_list))
        for i in range(len(R0_list)):
            dE_list[i] = dE_per_monomer(a, R_cell, R0_list[i], ks, ktheta, kl) * (n)
            
        popt, pcov = curve_fit(linear, x_list, dE_list)
        m = popt[0]
        l = m / (r_mono * (n))
        l_err = np.sqrt(np.diag(pcov))[0] / (r_mono * (n))
        
        theoryline = linear(x_fit, *popt)
        plt.plot(x_fit, theoryline, color=thiscolor, linestyle=':', lw = 1, label=r'$k_s = {:.2f}\,k_BT,\, k_{{\theta}} = {:.2f}\,k_BT$, $l_p = {:.2f} \pm {:.2f}\,\mathrm{{nm}}$ (theory)'.format(ks, ktheta, l, l_err), markersize = 3)
        
        print("lp (theory) = {:.2f} +/- {:.2f} nm".format(l, l_err))
        
        print("-" * 10)
    
# plt.xlabel(r'$(1/R_{\rm cell} - 1/R_0)^2\,(\mathrm{nm}^{-2})$')
plt.xlabel(r'$\left(\frac{1}{R_{\rm cell}} - \frac{1}{R_0}\right)^2\,(\mathrm{nm}^{-2})$', fontsize=14)
plt.ylabel(r'$\Delta E_{\mathrm{bend}}\,(k_BT)$', fontsize=14)

# plt.xscale('log')
# plt.yscale('log')

plt.legend(fontsize=8, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

plt.savefig("compare_Ebend_vs_R0.pdf")
