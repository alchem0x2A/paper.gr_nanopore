from utils import *
import numpy
import matplotlib.pyplot as plt
import os, os.path
from scipy.constants import k, e, electron_volt, epsilon_0


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

pairs = [("c_p", 1), ("c_n", -1)]
T = 300


plt.style.use("science")

for conc in [0.001]:
    file_name = os.path.join(out_path, file_template.format(conc))
    data = get_data(file_name)  # already in nm
    X = data[:, 0].reshape(*pixels); Y = data[:, 1].reshape(*pixels)
    x = data[:, 0]; y = data[:, 1]
    # cond = numpy.where(x)
    # print(X[1, 1], Y[1, 1])
    for quant, z in pairs:
        for Vg in Vg_plot:
            print(get_col_index(Vg, quant))
            c = data[:, get_col_index(Vg, quant)]
            c[numpy.isnan(c)] = 0
            print(numpy.max(c), numpy.min(c))
            v = numpy.nan_to_num(data[:, get_col_index(Vg, "V")])
            v0 = numpy.mean(v[y>19.0])
            mu = (k * T * numpy.log(c / (conc * ratio_conc)) + z * e * v) / electron_volt
            mu = (z * e * v) / electron_volt
            # mu = (k * T * numpy.log(c / (conc * ratio_conc))) / electron_volt
            mu0 = numpy.mean(mu[y>19.5])
            mu = mu - mu0
            D = mu.reshape(*pixels)
            D[numpy.isinf(D)] = 0
            print(numpy.max(D), numpy.min(D))
            plt.cla()
            fig = plt.figure(figsize=(2.8, 2.8))
            ax = fig.add_subplot(111)
            # if z > 0:
                # vmin = -0.10; vmax = 0.01
            # else:
                # vmin = 0.008; vmax = 0.011
            mesh = plot_data(ax, X, Y, D)
            mesh.set_cmap("jet")
            ax.set_title("{0} mol/L-{1} V-{2}".format(conc,
                                                            Vg,
                                                            quant
            ))
            ax.set_xlabel("$r$ (nm)")
            ax.set_ylabel("$z$ (nm)")
            add_graphene(ax, R_p=10)
            fig.colorbar(mesh, fraction=0.03)
            fig.tight_layout()
            outfile = os.path.join(plot_path,
                                     "mu-{0}-{1}-{2}.svg".format(conc,
                                                              Vg,
                                                              quant))
            print(outfile)
            fig.savefig(outfile)
            
    
    
