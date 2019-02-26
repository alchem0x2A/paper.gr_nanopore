from utils import *
import matplotlib
import numpy
import matplotlib.pyplot as plt
import os, os.path
from os.path import dirname, join, exists
from scipy.constants import pi, hbar, e
from scipy.interpolate import griddata
vf = 1.1e6

curdir = dirname(__file__)

"""
Plot Vg surface contour
"""

# Use C/m^2
def delta_phi_gr(sigma):
    fac = hbar * vf / e * numpy.sqrt(pi * numpy.abs(sigma) / e)
    return fac * numpy.sign(sigma)

quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
units = ("V", "mol/m$^{3}$", "mol/m$^{3}$", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")
# get the index of column in the new matrix
def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 1 + len_quant * idx_V + idx_quant

Vg_true = []
Vg_all = [0.001, *numpy.arange(0.025, 0.251, 0.025)]

file_template = "{0}.npy"
sigma_file_template = "sigma_{0}.txt"

concentrations = (0.0001, 0.0002, 0.0005,
                  0.001, 0.002, 0.005,
                  0.01, 0.02, 0.05, 0.1)

ratio_conc = 10**3

out_path = join(curdir, "../data/FEM/concentration/1D/")
plot_path = join(curdir, "../img/")

res = []
for i, conc in enumerate(concentrations):
    sigma_file = join(out_path, sigma_file_template.format(conc))
    sigma_data = numpy.genfromtxt(sigma_file, comments="%")
    V_ = sigma_data[:, 0]; sig_ = -sigma_data[:, 1]
    Vg_ = V_ + delta_phi_gr(sig_)
    [res.append((Vg_[i], conc, V_[i])) for i in range(len(V_))]

res = numpy.array(res)

v = res[:, 0]
lambda_d = Debye_length(res[:, 1]) / 1e-9
psi_g = res[:, 2]

res[:, 1] = lambda_d

v_uni = numpy.linspace(0.0, 1.5, 128)
l_uni = numpy.linspace(min(lambda_d), max(lambda_d), 128)
vv, ll = numpy.meshgrid(v_uni, l_uni)
z_uni = griddata(res[:, :2], psi_g,
                 (vv, ll), method="cubic")

r_p = 20

# print(res)

plt.figure(figsize=(2.8, 2.3))
plt.style.use("science")
# plt.yscale("log")

plt.xlabel("$V_{\mathrm{G}}$ (V)")
plt.ylabel("$c_{0}$ (mol$\cdot$L$^{-1}$)")
plt.scatter(v[psi_g>0.001], lambda_d[psi_g>0.001] / r_p, c=psi_g[psi_g>0.001],
            s=10,
            cmap="rainbow", alpha=0.23)
cs = plt.contour(vv, ll / r_p, z_uni, cmap="rainbow",
                 levels=Vg_all, vmin=0.0, vmax=0.25)
norm = matplotlib.colors.Normalize(vmin=0.0, vmax=0.25)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cs.cmap)
sm.set_array([])
plt.colorbar(mappable=sm)

plt.xlim(0, 1.25)



'''
for i, vg in enumerate(Vg_all[1: ]):
    cond = numpy.where(numpy.abs(res  - vg) < 1e-3)[0]  # row numbers
    plt.plot(res[cond, 0], res[cond, 1], "o-",
             label=vg,
             color=cmap[i + 1])
''' 
# for v in Vg_all:
    # plt.axvline(x=v)

# plt.xlim(0, 1.25)
plt.legend()
plt.tight_layout()
plt.savefig(join(plot_path, "countour_VG.svg"))
