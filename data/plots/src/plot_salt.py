from utils import *
import numpy
import matplotlib.pyplot as plt
import os, os.path
from scipy.constants import pi, hbar, e, N_A, k, epsilon_0
from scipy.interpolate import interp1d

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
sigma_file_template = "sigma_{0}_out.txt"

salts = ["NaCl", "LiCl", "KCl", "CaCl2",
         "MgSO4", "K2SO4", "K3FeCN"]

rates = [0.95, 1.0, 0.9, 0.35, 0.6, 0.55, 0.1]

ratio_conc = 10**3

# out_path = "../result/salt/"
out_path = "../data/FEM/salt/20/"
plot_path = "../img/"
plt.style.use("science")

# get the index of column in the new matrix
def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 1 + len_quant * idx_V + idx_quant


avg_tflux = []

quant_flux = ("zflux_cp", "zflux_cn")

fig = plt.figure(figsize=(4, 4))
for salt in salts:
    tflux = []
    file_name = os.path.join(out_path, file_template.format(salt))
    print(file_name)
    data = numpy.load(file_name)
    data[:, 0] /= 1e-9
    r = data[:, 0]
    for V in Vg_all:
        avg = 0
        idx = get_col_index(V, "zflux_cp")
        avg += numpy.abs(numpy.mean(data[:, idx]))
        idx = get_col_index(V, "zflux_cn")
        avg += numpy.abs(numpy.mean(data[:, idx]))
        tflux.append(avg)
    avg_tflux.append(1 - tflux / tflux[0])  # averaged
    sigma_file = os.path.join(out_path, sigma_file_template.format(salt))
    sigma_data = numpy.genfromtxt(sigma_file, comments="%")
    V_ = sigma_data[:, 0]; sig_ = -sigma_data[:, 1]
    Vg_true.append(V_ + delta_phi_gr(sig_))
    
Vg_true = numpy.array(Vg_true)
avg_tflux = numpy.array(avg_tflux)

fig = plt.figure(figsize=(2.5, 2.5))
plt.style.use("science")

for i, salt in enumerate(salts):
    # if salt == "CaCl2":
        # rate = 0.4
    # else:
        # rate = 1.0
    rate =  rates[i]
    x = Vg_true[i, :]
    y = avg_tflux[i, :]
    print(salt, y[x<=1.25][-1] * rate)
    xx = numpy.linspace(x.min(), 1.25, 128)
    yy = interp1d(x, y, kind="cubic")(xx)
    numpy.savetxt(fname="{}-FEM-curve.csv".format(salt),
                  X=numpy.vstack([xx, yy * rate]).T,
                  header="V_G (V), xi",
                  delimiter=",",
                  comments="#")
    plt.plot(xx, yy * rate, "-",
             markersize=5,
             label="{}".format(salt))

# plt.yscale("log")
plt.xlabel("$V_{g}$ (V)")
plt.xlim(0, 1.25)
# plt.ylim(-0.2, 1.0)
plt.ylabel("Rectification")
plt.legend()

plt.savefig(os.path.join(plot_path, "tflux_salt.svg"))


d_debye = Debye_length(numpy.array(concentrations)) / 1e-9
rect = {"KCl": ( 0.268, 0.055, 0.623),
        "NaCl": (0.563, 0.110, 0.703),
        "LiCl": (0.459, 0.096, 0.73),
        "CaCl2": (0.214, 0.049, 0.268),
        "K2SO4": (0.153, 0.052, 0.089),
        "MgSO4": (0.285, 0.069, 0.207),
        "K3FeCN": (-0.036, 0.040, -0.056)}

# rect = 1 - avg_tflux[:, 6] / avg_tflux[:, 0]
# salts = ["KCl"] + salts
# print(salts)
# rect = [0.82] + list(rect)


plt.figure(figsize=(2.5, 2.5))
plt.bar(numpy.arange(len(salts))-0.125, height=[rect[k][0] for k in salts],
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
print(LD)

from scipy.stats import linregress
regres = linregress(LD, [rect[s][-1] for s in salts])
k = regres.slope; b = regres.intercept; r = regres.rvalue
print(k, b, r)

plt.figure(figsize=(2.5, 2.5))
# plt.plot(LD, [rect[s][0] for s in salts], "o", label="Exp", alpha=0.8)
plt.errorbar(LD, [rect[s][0] for s in salts],  yerr=[rect[s][1] for s in salts], fmt="o", alpha=0.8, label="Exp")
plt.plot(LD, [rect[s][-1] for s in salts], "D", label="FEM", alpha=0.8)
xx = numpy.linspace(10, 35)
yy = k * xx  + b
plt.plot(xx, yy, "--")
# plt.xscale("log")
plt.xlabel("LD")
plt.xlim(12, 33)
plt.ylim(-0.08, 0.8)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(plot_path, "rect_comparison.svg"))

