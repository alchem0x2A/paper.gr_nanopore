import numpy
import matplotlib.pyplot as plt
import os
from os.path import dirname, join, exists
curdir = dirname(__file__)

"""
Plot experimental rectification for all salts
"""

file_templ = join(curdir, "../data/exp/{}.csv")
def read(name):
    file_name = file_templ.format(name)
    err = None
    with open(file_name) as f:
        err = float(f.readline().strip())
    data = numpy.genfromtxt(file_name,
                            skip_header=1,
                            delimiter=",")
    data[:, 1] -= data[data[:, 0] == 0, 1]
    return data, err

def convert(eta, delta):
    assert delta > 0
    xi_inv = (1 - eta) / (delta * eta + 1)
    return 1 - xi_inv

if __name__ == "__main__":
    fig = plt.style.use("science")
    plt.figure(figsize=(2.8, 2.8), facecolor="w")
    plt.axvline(x=0, ls="--", color="#00ffee")
    plt.axhline(y=0, ls="--", color="#00ffee")

    delta0 = 2.2
    for name in ["KCl", "NaCl", "LiCl",
                 "MgSO4", "K2SO4", "CaCl2",
                 "K3FeCN"]:
        data, err = read(name)
        v = data[:, 0]
        eta = data[:, 1]
        xi = convert(eta, delta=delta0)
        err_lo = xi - convert(eta - err / 2, delta=delta0)
        err_hi = convert(eta + err / 2, delta=delta0) - xi
        plt.errorbar(data[:, 0], xi,
                     yerr=(err_lo / 2, err_hi / 2),
                     # yerr=err / 2,
                     fmt="o-", alpha=0.5, label=name,
                     markersize=4)
        # plt.fill_between(data[:, 0], data[:, 1] - err / 2,
                         # data[:, 1] + err / 2,
                         # fmt="o",
                         # alpha=0.5, label=name)
    plt.ylim(-0.1, 0.9)
    plt.xlim(-1.25, 1.25)
    plt.xlabel("$V_{\\mathrm{G}}$ (V)")
    plt.ylabel("Rectification")
    plt.legend()
    plt.savefig(join(curdir, "../img/rect_all.svg".format(name)))
