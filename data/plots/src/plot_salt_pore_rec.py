from utils import *
import numpy
import matplotlib.pyplot as plt
import os
from os.path import dirname, join, exists
from scipy.constants import pi, hbar, e
import numpy
from scipy.interpolate import griddata, interp2d, interp1d
import pickle

"""
Plot rectification for all salts based on FEM simulation
"""

curdir = dirname(__file__)

vf = 1.1e6

salts =  ["NaCl", "LiCl", "KCl", "CaCl2",
         "MgSO4", "K2SO4", "K3FeCN"]

# rates = [0.95, 0.96, 0.92, 0.35, 0.6, 0.55, 0.1]

# Use C/m^2
def delta_phi_gr(sigma):
    fac = hbar * vf / e * numpy.sqrt(pi * numpy.abs(sigma) / e)
    return fac * numpy.sign(sigma)

def get_col_index(Vg, quantity):
    idx_V = Vg_all.index(Vg)
    idx_quant = quantities.index(quantity)
    len_quant = len(quantities)
    return 1 + len_quant * idx_V + idx_quant

r_list = list(numpy.arange(7.5, 35.1, 2.5))

def get_xi_salt(salt):
    assert salt in ("KCl", "NaCl", "LiCl",
                    "MgSO4", "K2SO4", "CaCl2", "K3FeCN")
    fname = join(curdir, "../data/FEM/salt/radius/{}_combined.csv").format(salt)
    if os.path.exists(fname):
        print(salt)
        res = numpy.genfromtxt(fname, delimiter=",", comments="%")
    else:
        Vg_all = [0.001, *numpy.arange(0.025, 0.301, 0.025)]
        Vg_true = []

        file_template = join(curdir,
                             "../data/FEM/salt/radius/{0}",
                             "{1}.npy")
        sigma_file_template = join(curdir,
                                   "../data/FEM/salt/radius/{0}",
                                   "sigma_{1}.txt")
        # get the index of column in the new matrix
        quantities = ("V", "c_p", "c_n", "zflux_cp", "zflux_cn")
        units = ("V", "mol/m$^{3}$", "mol/m$^{3}$", "mol/(m$^{2}$*s)", "mol/(m$^{2}$*s)")
        avg_tflux = []
        quant_flux = ("zflux_cp", "zflux_cn")
        res = []
        for r_g in r_list:
            tflux = []
            # conc_cs = str(conc).split(".")[-1]
            # file_name = join(out_path, file_template.format(conc))
            file_name = file_template.format(salt, int(r_g))
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
            sigma_file = sigma_file_template.format(salt, int(r_g))
            # sigma_file = join(out_path, sigma_file_template.format(conc))
            sigma_data = numpy.genfromtxt(sigma_file, comments="%")
            V_ = sigma_data[:, 0]; sig_ = -sigma_data[:, 1]
            Vg_ = V_ + delta_phi_gr(sig_)
            [res.append((Vg_[i], r_g - 5, rec_[i])) for i in range(len(V_)) if Vg_[i] < 1.5]
        res = numpy.array(res)
        numpy.savetxt(fname=fname,
                      X=res, delimiter=",", header="Vg,r,xi")
    # print(res)
    return res

def get_1D(salt,
           r_list=numpy.arange(2.5, 30.1, 2.5),
           interpolate=True):
    res = get_xi_salt(salt)
    v = res[:, 0]
    r = res[:, 1]
    xi = res[:, 2]
    y_interp = []
    x_raw = []
    y_raw = []
    if not interpolate:
        raise NotImplementedError
    # interpolate
    x = numpy.linspace(0, 1.25, 15)
    for r_ in r_list:
        y = interp1d(v[r==r_], xi[r==r_], fill_value="extrapolate")(x)
        y_interp.append(y)
        x_raw.append(v[r==r_])
        y_raw.append(xi[r==r_])
    return x, y_interp, x_raw, y_raw

r_list =  numpy.arange(2.5, 30.1, 2.5)
cmap_name = "rainbow_r"
colors = plt.cm.get_cmap("rainbow_r")(numpy.linspace(0, 1, len(r_list)))

def plot_salt_trend(salt):
    fig = plt.figure(figsize=(2.8, 2.3))
    plt.style.use("science")
    x, y, x_raw, y_raw = get_1D(salt)
    # print(x, y, x_raw, y_raw)
    # ratio = rates[salts.index(salt)]
    for i, r_ in enumerate(numpy.arange(2.5, 30.1, 2.5)):
        plt.plot(x, y[i], color=colors[i])
        plt.plot(x_raw[i], y_raw[i], "o",
                 markersize=5,
                 color=colors[i])
    plt.scatter([0, 0], [0, 0], s=0,
                c=[2.5, 35], cmap="rainbow_r")
    plt.axhline(y=0, ls="--", color="grey")
    cb = plt.colorbar(shrink=0.6)
    cb.ax.set_title("$r_{\\mathrm{G}}$ (nm)")
    plt.ylim(-0.10, 1)
    plt.xlim(0, 1.25)
    plt.savefig(join(curdir, "../img/test-salt-rect-{}.svg").format(salt))
    
def plot_rec_pore_dist():
    file_pore_dist = join(curdir, "../data/exp/pore-dist.csv")
    exp_data = numpy.genfromtxt(file_pore_dist,
                                delimiter=",")
    rg = exp_data[:, 0] / 2 + 2.5
    w = exp_data[:, 1]
    w = w / numpy.sum(w)
    fig = plt.figure(figsize=(2.8, 2.3))
    plt.style.use("science")
    results = []
    for salt in salts:
        x, y, *_ = get_1D(salt)
        print(x)
        y = numpy.array(y); print(y.shape)
        s = 0
        for r_, w_ in zip(rg, w):
            if r_ in r_list:
                y_ = y[r_list==r_, :]
                s += y_ * w_
        s = s[0]
        xx = numpy.linspace(0, 1.25, 1024)
        ss = interp1d(x, s, kind="cubic")(xx)
        # plt.plot(x, ratio * s, label=salt)
        plt.plot(xx, ss, label=salt)
        print(ss)
        results.append("{0}, {1}\n".format(salt, ss[xx==1.25]))
    plt.legend()
    plt.ylim(-0.10, 1)
    plt.xlim(0, 1.25)
    plt.xlabel("$V_{\\mathrm{G}}$ (V)")
    plt.ylabel("$\\xi$")
    plt.savefig(join(curdir, "../img/FEM-salts-pore-dist.svg"))
    with open(join(curdir,
                   "../data/FEM/salt/radius/",
                   "final-ratio.csv"), "w") as f:
        f.writelines(results)
        

if __name__ == "__main__":
    # Plot salts
    for salt in ("KCl", "NaCl", "LiCl",
                 "MgSO4", "K2SO4", "CaCl2", "K3FeCN"):
        plot_salt_trend(salt)

    # Combine exp values
    plot_rec_pore_dist()




'''

data_pore = numpy.genfromtxt(join(curdir, "../data/exp/pore-dist.csv"),
                             delimiter=",")
r_exp = data_pore[:, 0] + 2.5
w_exp = data_pore[:, 1]
w_exp = w_exp / numpy.sum(w_exp)            # Frequencies

r_g = 20
conc = 1.1e-4
lambda_d = Debye_length(conc) / 1e-9
print(lambda_d)

fig = plt.figure(figsize=(2.8, 2.3))
plt.style.use("science")

v = numpy.linspace(0, 1.25, 128)
# l = numpy.ones_like(v) *
l = lambda_d / r_g
xi_sim = func_rec_interp(v, l)

conc = 2e-4
lambda_d = Debye_length(conc) / 1e-9
l_exp = lambda_d / r_exp
xi_exp = func_rec_interp(v, l_exp)
xi_exp = numpy.dot(w_exp, xi_exp)




plt.plot(v, xi_sim.flat, label="Single pore")
plt.plot(v, xi_exp.flat, label="With pore distribution")
plt.xlabel("$V_{\\mathrm{G}}$ (V)")
plt.ylabel("$\\xi$")
plt.legend()

plt.savefig(join(curdir, "../img/simple-rect-pore-dist.svg"))
'''
