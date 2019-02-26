from utils import *
import numpy
import matplotlib.pyplot as plt
import os, os.path


pixels = (512, 512)

quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
units = ("V", "mol/L", "mol/L", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")

Vg_all = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2] + list(numpy.arange(0.25, 1.35, 0.1))

file_template = "{0}.npy"

concentrations = (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)

ratio_conc = 10**3

# get the index of column in the new matrix
def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 2 + len_quant * idx_V + idx_quant


out_path = "../result/concentration/2D/"
plot_path = "../plot/concentration"

# Vg_plot = (0.001, 0.15)
Vg_plot = (0.15,)

plt.style.use("science")

for conc in [0.001]:
    file_name = os.path.join(out_path, file_template.format(conc))
    data = get_data(file_name)  # already in nm
    X = data[:, 0].reshape(*pixels); Y = data[:, 1].reshape(*pixels)
    print(X[1, 1], Y[1, 1])
    for quant in quantities:
        vmins = []; vmaxs = []
        for Vg in Vg_plot:
            print(get_col_index(Vg, quant))
            D = data[:, get_col_index(Vg, quant)]
            vmins.append(numpy.min(D)); vmaxs.append(numpy.max(D))
        for Vg in Vg_plot:
            D = data[:, get_col_index(Vg, quant)].reshape(*pixels)
            plt.cla()
            fig = plt.figure(figsize=(2.8, 2.8))
            ax = fig.add_subplot(111)
            if quant in ("zflux_cp", "zflux_cn"):
                vmax = 0
            else:
                vmax = None
            mesh = plot_data(ax, X, Y, D,
                             vmin=-0.04,
                             vmax=0)
            if quant in ("zflux_cp", "zflux_cn"):
                mesh.set_cmap("jet_r")
            ax.set_title("{0} mol/L-{1} V-{2} ({3})".format(conc,
                                                            Vg,
                                                            quant,
                                                            units[quantities.index(quant)]
            ))
            ax.set_xlabel("$r$ (nm)")
            ax.set_ylabel("$z$ (nm)")
            add_graphene(ax, R_p=10)
            fig.colorbar(mesh, fraction=0.03)
            fig.tight_layout()
            outfile = os.path.join(plot_path,
                                     "{0}-{1}-{2}.svg".format(conc,
                                                              Vg,
                                                              quant))
            print(outfile)
            fig.savefig(outfile)
            
    
    
