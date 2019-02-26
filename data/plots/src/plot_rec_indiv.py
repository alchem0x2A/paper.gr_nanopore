import numpy
import matplotlib.pyplot as plt
import os
from os.path import dirname, join, exists
curdir = dirname(__file__)

"""
Plot experimental rectification of individual salts
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

# def convert(eta, delta):
    # assert delta > 0
    # xi_inv = (1 - eta) / (delta * eta + 1)
    # return 1 - xi_inv

if __name__ == "__main__":
    fig = plt.style.use("science")

    for name in ["KCl", "NaCl", "LiCl",
                 "MgSO4", "K2SO4", "CaCl2",
                 "K3FeCN"]:
        plt.cla()
        plt.figure(figsize=(2, 2), facecolor="w")
        plt.axvline(x=0, ls="--", color="grey")
        plt.axhline(y=0, ls="--", color="grey")
        data, err = read(name)
        v = data[:, 0]
        eta = data[:, 1]
        plt.errorbar(data[:, 0], eta,
                     yerr=err / 2,
                     fmt="o-", alpha=0.5, label=name,
                     markersize=4)
        plt.ylim(-0.1, 0.7)
        plt.xlim(-1.25, 1.25)
        plt.xlabel("$V_{\\mathrm{G}}$ (V)")
        plt.ylabel("$\\eta$")
        plt.title(name)
        plt.savefig(join(curdir, "../img/rect_{}.svg".format(name)))
