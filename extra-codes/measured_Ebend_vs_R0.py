import numpy as np
import matplotlib.pyplot as plt
import rcparams
from scipy.optimize import curve_fit

colorlist = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'cyan', 'magenta', 'olive', 'lime', 'teal', 'lavender', 'salmon', 'turquoise', 'tan', 'gold', 'orchid', 'skyblue', 'lightgreen', 'lightblue', 'lightcoral', 'lightpink', 'lightgray', 'lightcyan', 'lightyellow', 'lightpurple', 'lightbrown', 'darkred', 'darkblue', 'darkgreen', 'darkorange', 'darkpurple', 'darkbrown', 'darkpink', 'darkgray', 'darkcyan', 'darkmagenta', 'darkolive', 'darklime', 'darkteal', 'darklavender', 'darksalmon', 'darkturquoise', 'darktan', 'darkgold', 'darkorchid', 'darkskyblue', 'darklightgreen', 'darklightblue', 'darklightcoral', 'darklightpink', 'darklightgray', 'darklightcyan', 'darklightyellow', 'darklightpurple', 'darklightbrown']


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


headers, datasets = read_data("measured_Ebend_vs_R0.dat")

Ebind_list, ktheta_list = headers

Ebind_banlist = []
Ebind_banlist = np.array(Ebind_banlist)

R_cell = 350
kBT0 = 310

r_mono = 2.5
n = 20

plt.figure(figsize=(5, 5), constrained_layout=True)

for Ebind_i in range(len(Ebind_list)):

    Ebind = Ebind_list[Ebind_i]
    ktheta = ktheta_list[Ebind_i]
    thiscolor = colorlist[Ebind_i % len(colorlist)]
    # thiscolor = plt.rcParams['axes.prop_cycle'].by_key()['color'][Ebind_i % len(plt.rcParams['axes.prop_cycle'].by_key()['color'])]

    if Ebind not in Ebind_banlist:

        Ebind = Ebind / kBT0
        ktheta = ktheta / kBT0

        print("Ebind = {:.2f} kBT, ktheta = {:.2f} kBT".format(Ebind, ktheta))

        R0_list, Eb_list = datasets[Ebind_i]

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

        plt.plot(x_list, Eb_list, 'o', color=thiscolor, markersize = 3, label=r'$E_{{\mathrm{{bind}}}} = {:.2f}\,k_BT,\, k_{{\theta}} = {:.2f}\,k_BT$'.format(Ebind, ktheta))
        plt.plot(x_fit, y_fit, color=thiscolor, linestyle='--', lw = 1, label=r'$l_p = {:.2f} \pm {:.2f}\,\mathrm{{nm}}$'.format(lp, lp_err))

        # print("Measured slope = {:.2f}".format(m))
        print("Persistence length = {:.2f} +/- {:.2f} nm".format(lp, lp_err))
    
# plt.xlabel(r'$(1/R_{\rm cell} - 1/R_0)^2\,(\mathrm{nm}^{-2})$')
plt.xlabel(r'$\left(\frac{1}{R_{\rm cell}} - \frac{1}{R_0}\right)^2\,(\mathrm{nm}^{-1})$', fontsize=14)
plt.ylabel(r'$E_{\mathrm{bend}}\,(k_BT)$', fontsize=14)

# plt.xscale('log')
# plt.yscale('log')

plt.legend(fontsize=8)

plt.savefig("measured_Ebend_vs_R0.pdf")
