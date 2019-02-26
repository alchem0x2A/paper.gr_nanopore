from utils import *
import numpy
import matplotlib.pyplot as plt
import os, os.path
from scipy.constants import pi, hbar, e
vf = 1.1e6

# Use C/m^2
def delta_phi_gr(sigma):
    fac = hbar * vf / e * numpy.sqrt(pi * numpy.abs(sigma) / e)
    return fac * numpy.sign(sigma)


quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
units = ("V", "mol/m$^{3}$", "mol/m$^{3}$", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")

# Vg_all = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2] + list(numpy.arange(0.25, 1.35, 0.1))
Vg_all = [0.001, *numpy.arange(0.025, 0.251, 0.025)]

Vg_true = []

file_template = "0.{0}.npy"
# file_template = "0_{0}_out_1Ddata.txt"
sigma_file_template = "sigma_0.{0}.txt"

concentrations = (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)
# concentrations = (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)
radii = (5, 10, 15, 20)

ratio_conc = 10 ** 3

out_path = "../result/radius/{0}/1D/"
plot_path = "../plot/radius/1D/"

# get the index of column in the new matrix
def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 1 + len_quant * idx_V + idx_quant

conc_plot = [0.001]

for rad in [10]:    
    avg_tflux = []
    avg_nflux = []
    avg_pflux = []

    quant_flux = ("zflux_cp", "zflux_cn")
    for conc in [0.001]:
        conc_base = str(conc).split(".")[-1]
        tflux = []; nflux = []; pflux = []
        file_name = os.path.join(out_path.format(rad), file_template.format(conc_base))
        # data = numpy.genfromtxt(file_name, comments="%")
        data = numpy.load(file_name)
        data[numpy.isnan(data)] = 0
        data[:, 0] /= 1e-9
        r = data[:, 0]
        for V in Vg_all:
            avg = 0
            idx = get_col_index(V, "zflux_cp")
            avg += numpy.abs(numpy.mean(data[:, idx]))
            pflux.append(numpy.abs(numpy.mean(data[:, idx])))
            idx = get_col_index(V, "zflux_cn")
            avg += numpy.abs(numpy.mean(data[:, idx]))
            nflux.append(numpy.abs(numpy.mean(data[:, idx])))
            # print(avg)
            tflux.append(avg)
        avg_tflux.append(tflux)
        avg_pflux.append(pflux)
        avg_nflux.append(nflux)
            # correction for Vg
        sigma_file = os.path.join(out_path.format(rad), sigma_file_template.format(conc_base))
        sigma_data = numpy.genfromtxt(sigma_file, comments="%")
        V_ = sigma_data[:, 0]; sig_ = -sigma_data[:, 1]
        Vg_true.append(V_ + delta_phi_gr(sig_))

    avg_tflux = numpy.array(avg_tflux)
    avg_nflux = numpy.array(avg_nflux)
    avg_pflux = numpy.array(avg_pflux)

    fig = plt.figure(figsize=(2.5, 2.5))
    plt.style.use("science")
    
    col = plt.cm.rainbow(numpy.linspace(0, 1, 7))

    for i, conc in enumerate(conc_plot):
        # print(avg_tflux[i, :])
        plt.plot(Vg_true[i], numpy.abs(avg_tflux[i, :]), "-o",
                 markersize=5,
                 label="{} mol/L$^3$".format(conc))
        plt.plot(Vg_true[i], numpy.abs(avg_nflux[i, :]), "-o",
                 markersize=5,
                 label="{} mol/L$^3$".format(conc), alpha=0.5)
        plt.plot(Vg_true[i], numpy.abs(avg_pflux[i, :]), "-o",
                 markersize=5,
                 label="{} mol/L$^3$".format(conc), alpha=0.5)

    # plt.yscale("log")
    plt.xlabel("True $V_{g}$ (V)")
    plt.xlim(0, 1.25)
    # plt.ylim(2e-4, 3e1)
    plt.ylabel("Average total flux (mol/(m$^{-2}$*s))")
    plt.legend()
    plt.title("{} nm".format(rad))

    plt.savefig(os.path.join(plot_path, "rect_Vg_{}.svg".format(rad)))


    
    
    
    
    
