from utils import *
import numpy
import matplotlib.pyplot as plt
import os, os.path

pixels = (512, 512)

quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
units = ("V", "mol/m$^{3}$", "mol/m$^{3}$", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")

Vg_all = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2] + list(numpy.arange(0.25, 1.35, 0.1))

file_template = "{0}.npy"

concentrations = (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)

ratio_conc = 10**3

out_path = "../result/concentration/1D/"
plot_path = "../plot/concentration/1D"

# get the index of column in the new matrix
def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 1 + len_quant * idx_V + idx_quant

Vg_plot = Vg_all[: 10]

plt.style.use("science")
fig = plt.figure(figsize=(4, 4))
for conc in concentrations:
    file_name = os.path.join(out_path, file_template.format(conc))
    data = numpy.load(file_name)
    data[:, 0] /= 1e-9
    r = data[:, 0]
    for quant in quantities:
        plt.cla()
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(111)
        for V in Vg_plot:
            y = data[:, get_col_index(V, quant)]
            plt.plot(r, y, label="{:.2f} V".format(V))
        ax.set_xlabel("$r$ (nm)")
        ax.set_ylabel("{0} ({1})".format(quant,
                                         units[quantities.index(quant)]))
        ax.set_title("c0 = {}, {}".format(conc, quant))
        ax.legend()
        fig.tight_layout()
        outfile = os.path.join(plot_path,
                               "1D-{}-{}.pdf".format(conc, quant))
        print(outfile)
        fig.savefig(outfile)
    
    
    
    
