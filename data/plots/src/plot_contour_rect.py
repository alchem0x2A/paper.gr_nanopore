from utils import *
import numpy
import matplotlib.pyplot as plt
import os
from os.path import dirname, join, exists
from scipy.constants import pi, hbar, e
import numpy
from scipy.interpolate import griddata

curdir = dirname(__file__)

"""
Plot contour plot of surface rectification
"""


vf = 1.1e6


# Use C/m^2
def delta_phi_gr(sigma):
    fac = hbar * vf / e * numpy.sqrt(pi * numpy.abs(sigma) / e)
    return fac * numpy.sign(sigma)

pixels = (512, 512)

quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
units = ("V", "mol/m$^{3}$", "mol/m$^{3}$", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")

# Vg_all = [0.001, 0.025, 0.05, 0.1, 0.15, 0.2] + list(numpy.arange(0.25, 1.35, 0.1))

Vg_all = [0.001, *numpy.arange(0.025, 0.251, 0.025)]
Vg_true = []

file_template = "{0}.npy"
sigma_file_template = "sigma_{0}.txt"

# concentrations = (0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)
concentrations = (0.0001, 0.0002, 0.0005,
                  0.001, 0.002, 0.005,
                  0.01, 0.02, 0.05, 0.1)

ratio_conc = 10**3

out_path = join(curdir, "../data/FEM/concentration/1D/")
plot_path = join(curdir, "../img/")


# get the index of column in the new matrix
def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 1 + len_quant * idx_V + idx_quant


avg_tflux = []

quant_flux = ("zflux_cp", "zflux_cn")

res = []
for conc in concentrations:
    tflux = []
    # conc_cs = str(conc).split(".")[-1]
    file_name = join(out_path, file_template.format(conc))
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
    rec_ = 1 - tflux / tflux[0]
    # correction for Vg
    sigma_file = join(out_path, sigma_file_template.format(conc))
    sigma_data = numpy.genfromtxt(sigma_file, comments="%")
    V_ = sigma_data[:, 0]; sig_ = -sigma_data[:, 1]
    Vg_ = V_ + delta_phi_gr(sig_)
    [res.append((Vg_[i], conc, rec_[i])) for i in range(len(V_))]

res = numpy.array(res)

v = res[:, 0]
lambda_d = Debye_length(res[:, 1]) / 1e-9
rec = res[:, 2]
cond = numpy.where(rec > 0)[0]
print(cond)

res[:, 1] = lambda_d

v_uni = numpy.linspace(0.0, 1.25, 128)
l_uni = numpy.linspace(min(lambda_d), 35, 128)
vv, ll = numpy.meshgrid(v_uni, l_uni)
z_uni = griddata(res[:, :2], rec,
                 (vv, ll), method="cubic", fill_value=0)
# print(res[:, :2])
# print(rec)
# print(z_uni)

r_p = 20

fig = plt.figure(figsize=(2.8, 2.3))
plt.style.use("science")

# plt.scatter(v, lambda_d / r_p,
            # c=rec, s=10, cmap="rainbow",
            # alpha=0.25)

cs = plt.pcolor(vv * 1.085, ll * 1.0855 / r_p, z_uni,
                rasterized=True,
                cmap="rainbow")
plt.colorbar()




# plt.yscale("log")
# for i, conc in enumerate(concentrations):
    # xx = Vg_true[i]
    # yy = numpy.ones_like(Vg_true[i]) * conc
    # cc = avg_tflux[i, :]
    # plt.scatter(xx, yy ** -0.5, c=cc, cmap="rainbow")

plt.xlim(0.1, 1.25)
plt.ylim(0.1, 1.4)

plt.xlabel("True $V_{g}$ (V)")
plt.ylabel("Average total flux (mol/(m$^{-2}$*s))")
plt.tight_layout()

plt.savefig(join(plot_path, "rect_Vg_contour.svg"))

# plt.cla()

# numpy.save
# cond = numpy.where(vv == 1.25)
# plt.plot(l_uni * 1.25 / r_p, z_uni[:, -3])
# print(z_uni[:, -3].max())
# plt.xlim(0.1, 1.4)
# plt.show()
