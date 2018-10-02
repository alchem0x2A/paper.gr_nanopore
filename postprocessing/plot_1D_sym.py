from utils import *
import numpy
import matplotlib.pyplot as plt
from scipy.constants import e, epsilon_0, pi, k, N_A, hbar
import os, os.path

pixels = (512, 512)

quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
units = ("V", "mol/m$^{3}$", "mol/m$^{3}$", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")

# Vg_all = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2] + list(numpy.arange(0.25, 1.35, 0.1))

Vg_all = [0.001, *numpy.arange(0.025, 0.251, 0.025)]
print(Vg_all)
file_template = "{0}.npy"
sigma_file_template = "sigma_{0}.txt"

concentrations = (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)

ratio_conc = 10**3

out_path = "../result/radius/10/1D/"
plot_path = "../plot/radius/1D/"

# get the index of column in the new matrix
def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 1 + len_quant * idx_V + idx_quant

Vg_plot = Vg_all[: 10]

vf = 1.1e6

# Use C/m^2
def delta_phi_gr(sigma):
    fac = hbar * vf / e * numpy.sqrt(pi * numpy.abs(sigma) / e)
    return fac * numpy.sign(sigma)

plot_quantities = ("zflux_cp", "zflux_cn")
plt.style.use("science")

r0 = 10
for conc in [0.001]:
    file_name = os.path.join(out_path, file_template.format(conc))
    data = numpy.load(file_name)
    data[:, 0] /= 1e-9          
    r = data[:, 0]              # r distance
    for quant in plot_quantities:
        plt.cla()
        fig = plt.figure(figsize=(2.8, 2.2))
        ax = fig.add_subplot(111)
        lines = []
        for V in Vg_all:
            y = data[:, get_col_index(V, quant)]
            lines.append(plt.plot(numpy.hstack([-r[::-1], r]) / r0,
                                  numpy.hstack([y[::-1], y]))[0])
        sigma_file = os.path.join(out_path, sigma_file_template.format(conc))
        sigma_data = numpy.genfromtxt(sigma_file, comments="%")
        Vg_true = sigma_data[:, 0] + delta_phi_gr(-sigma_data[:, 1])  # True Vg
        print(Vg_true)
        for vg, line in zip(Vg_true, lines):
            line.set_color(plt.cm.rainbow(vg / 1.25))
        import matplotlib
        norm = matplotlib.colors.Normalize(vmin=0.0, vmax=1.25)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=plt.cm.rainbow)
        sm.set_array([])
        plt.colorbar(mappable=sm)
        # ax.set_xlabel("$r$ (nm)")
        # ax.set_ylabel("{0} ({1})".format(quant,
                                         # units[quantities.index(quant)]))
        # ax.set_title("c0 = {}, {}".format(conc, quant))
        fig.tight_layout()
        outfile = os.path.join(plot_path,
                               "1D-{}-{}.svg".format(conc, quant))
        print(outfile)
        fig.savefig(outfile)

    
    
    
    
