from utils import *
import numpy
import matplotlib.pyplot as plt
import os, os.path

pixels = (512, 512)

quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
units = ("V", "mol/m$^{3}$", "mol/m$^{3}$", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")

# Vg_all = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2] + list(numpy.arange(0.25, 1.35, 0.1))
Vg_all = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65]

file_template = "{0}.npy"

salts = ["NaCl", "LiCl",
         "MgSO4", "K2SO4", "K3FeCN"]

ratio_conc = 10**3

out_path = "../result/salt/"
plot_path = "../plot/salt"

# get the index of column in the new matrix
def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 1 + len_quant * idx_V + idx_quant


avg_tflux = []

quant_flux = ("zflux_cp", "zflux_cn")

plt.style.use("science")
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
    avg_tflux.append(tflux)

avg_tflux = numpy.array(avg_tflux)

fig = plt.figure(figsize=(4, 4))
plt.style.use("science")

for i, salt in enumerate(salts):
    plt.plot(Vg_all, avg_tflux[i, :], "-o",
             markersize=5,
             label="{}".format(salt))

plt.yscale("log")
plt.xlabel("$V_{g}$ (V)")
plt.ylabel("Average total flux (mol/(m$^{-2}$*s))")
plt.legend()

plt.savefig(os.path.join(plot_path, "tflux_salt.pdf"))


# d_debye = Debye_length(numpy.array(concentrations)) / 1e-9
rect = 1 - avg_tflux[:, 6] / avg_tflux[:, 0]
salts = ["KCl"] + salts
print(salts)
rect = [0.82] + list(rect)

plt.figure(figsize=(4, 4))
plt.bar(range(len(salts)), rect, 0.5)
plt.gca().set_xticks(range(len(salts)))
plt.gca().set_xticklabels(salts)
plt.ylabel("Rectification ratio")
plt.savefig(os.path.join(plot_path, "rect_salts.pdf"))

    
    
    
    
    
