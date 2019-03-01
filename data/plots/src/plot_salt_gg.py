from utils import *
import numpy
import matplotlib.pyplot as plt
import os, os.path
from os.path import join, dirname, exists
from scipy.constants import pi, hbar, e, N_A, k, epsilon_0

curdir = dirname(__file__)

vf = 1.1e6

# Use C/m^2
def delta_phi_gr(sigma):
    fac = hbar * vf / e * numpy.sqrt(pi * numpy.abs(sigma) / e)
    return fac * numpy.sign(sigma)

pixels = (512, 512)

quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
units = ("V", "mol/m$^{3}$", "mol/m$^{3}$", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")

# Vg_all = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2] + list(numpy.arange(0.25, 0.651, 0.1))
Vg_all = [0.001] + list(numpy.arange(0.025, 0.301, 0.025))
Vg_true = []

file_template = "{0}.npy"
sigma_file_template = "sigma_{0}.txt"

salts = ["NaCl", "LiCl", "KCl", "CaCl2",
         "MgSO4", "K2SO4", "K3FeCN"]

ratio_conc = 10**3

# out_path = "../result/salt/"
out_path = join(curdir, "../data/FEM/salt/20/")
plot_path = join(curdir, "../img")
plt.style.use("science")

d_debye = Debye_length(numpy.array(concentrations)) / 1e-9
rect = {"KCl": ( 0.268, 0.055, 0.832),
        "NaCl": (0.563, 0.110, 0.865),
        "LiCl": (0.459, 0.096, 0.885),
        "CaCl2": (0.214, 0.049, 0.216),
        "K2SO4": (0.153, 0.052, 0.284),
        "MgSO4": (0.285, 0.069, 0.430),
        "K3FeCN": (-0.016, 0.040, 0.027)}

# rect = 1 - avg_tflux[:, 6] / avg_tflux[:, 0]
# salts = ["KCl"] + salts
# print(salts)
# rect = [0.82] + list(rect)
def convert(delta, eta):
    assert delta > 0
    xi_inv = (1 - eta) / (delta * eta + 1)
    return 1 - xi_inv

plt.figure(figsize=(2.5, 2.5))
delta = 2.2
plt.bar(numpy.arange(len(salts))-0.125,
        height=[convert(delta, rect[k][0]) for k in salts],
        yerr=numpy.array(numpy.array([rect[k][1] for k in salts])),
        width=0.25, alpha=0.5, label="Exp")
plt.bar(numpy.arange(len(salts))+0.125, height=[rect[k][-1] for k in salts],
        width=0.25, alpha=0.5, label="FEM")
plt.gca().set_xticks(range(len(salts)))
plt.gca().set_xticklabels(salts, rotation=-45, ha="left")
plt.ylabel("Rectification")
plt.legend()
plt.savefig(os.path.join(plot_path, "rect_salts.svg"))

salts = ["NaCl", "LiCl", "KCl", "CaCl2",
         "MgSO4", "K2SO4", "K3FeCN"]
c_p = numpy.array([1, 1, 1, 1, 1, 2, 3]) * 1e-4 * ratio_conc
c_n = numpy.array([1, 1, 1, 2, 1, 1, 1]) * 1e-4 * ratio_conc
z_p = numpy.array([1, 1, 1, 2, 2, 1, 1])
z_n = numpy.array([1, 1, 1, 1, 2, 2, 3])

I = (z_p ** 2 * c_p + z_n ** 2 * c_n) / 2.0
eps = 78;
LD = (eps * epsilon_0 * k * 300 / (2 * N_A * e ** 2 * I)) ** 0.5 / 1e-9

from scipy.stats import linregress
regres = linregress(LD, [rect[s][1] for s in salts])
k = regres.slope; b = regres.intercept; r = regres.rvalue
print(k, b, r)

plt.figure(figsize=(2.5, 2.5))
# plt.plot(LD, [rect[s][0] for s in salts], "o", label="Exp", alpha=0.8)
plt.errorbar(LD, [convert(delta, rect[s][0]) for s in salts],  yerr=[rect[s][1] for s in salts], fmt="o", alpha=0.8, label="Exp")
plt.plot(LD, [rect[s][-1] for s in salts], "D", label="FEM", alpha=0.8)
# xx = numpy.linspace(3, 35)
# yy = k * xx  + b
# plt.plot(xx, yy, "--")
# plt.xscale("log")
plt.xlabel("LD")
plt.xlim(10, 33)
plt.ylim(-0.1, 1.0)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(plot_path, "rect_comparison_gg.svg"))

