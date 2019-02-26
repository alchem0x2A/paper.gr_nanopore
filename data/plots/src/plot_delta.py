import numpy
import matplotlib.pyplot as plt
import os
from os.path import dirname, join, exists
curdir = dirname(__file__)

"""
Plot estimated delta from experimental data
"""

salts = ["KCl", "NaCl", "LiCl", "CaCl2",
         "MgSO4", "K3FeCN6", "K2SO4"]
data = numpy.genfromtxt(join(curdir, "../data/exp/diffuse-pcte-salt-delta.csv"),
                        skip_header=2,
                        delimiter=",")
exclude = ["LiCl", "MgSO4"]

cond = [i for i, s in enumerate(salts) if s not in exclude]
cond_salts = [s for i, s in enumerate(salts) if s not in exclude]
c = data[cond, -3]
low = data[cond, -2]
up = data[cond, -1]
delta = 2.20

plt.figure(figsize=(3, 3 * 0.8))
plt.style.use("science")
plt.errorbar(range(len(cond)), y=c, yerr=[c - low, up - c],
             fmt="s", label="Experiment")
plt.axhline(y=delta, ls="--", label="Model")
# plt.plot(data[:, 1][cond] / 3600, data[:, 2][cond] / 1e-3)
# plt.plot(data[:, 1][data[:, 1] > 7200] / 3600, data[:, 2][data[:, 1] > 7200] / 1e-3)
plt.xticks(range(len(cond)))
plt.legend(loc=0)
plt.gca().set_xticklabels(cond_salts, rotation=-45, ha="left")
plt.xlim(-0.5, 5.5)
# plt.xlabel("t (h)")
plt.ylabel("$\\delta$")
plt.ylim(0, 10)
plt.savefig(join(curdir, "../img/delta-estimate.svg"))
