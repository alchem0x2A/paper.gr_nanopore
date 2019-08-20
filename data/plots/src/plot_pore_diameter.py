import numpy
import matplotlib.pyplot as plt
import os
from os.path import dirname, join, exists
from utils import Debye_length
import pickle
curdir = dirname(__file__)

"""
Plot experimental pore distribution
"""

file_name = join(curdir, "../data/exp/pore-dist.csv")
data = numpy.genfromtxt(file_name, delimiter=",")

fig = plt.style.use("science")
plt.figure(figsize=(2.5, 2.5), facecolor="w")
repeat = [d for d, n in data for _ in range(int(n))]
plt.hist(repeat, numpy.arange(5,66,5))
# plt.ylim(0, 0.04)
 # plt.yticks([0, 0.02, 0.04])
plt.xlabel("Pore Radius (nm)")
plt.ylabel("Frequencies")
plt.savefig(join(curdir, "../img/radius_dist.svg"))



"""
Exp pore distr with lambda_D
"""
lambda_D = Debye_length(c0=1.1e-4)
print(lambda_D)
rr = numpy.arange(10, 110)
xx = lambda_D / ((rr + 10) * 1e-9) * 1.75

fig = plt.style.use("science")

func_rec_interp = pickle.load(open(join(curdir,
                                   "../data/FEM/concentration/1D",
                                   "rect_2d_intep.pickle"), "rb"))

xi = func_rec_interp(1.25, xx)

plt.figure(figsize=(2.8, 2.3), facecolor="w")
repeat = [d for d, n in data for _ in range(int(n * (d /2) ** 2))]
plt.hist(repeat, numpy.arange(5,66,5), normed=True)
ax = plt.gca()
ax2 = ax.twiny()
ax3 = ax.twinx()
ax2.set_xticks(lambda_D / numpy.arange(0.5, 3.5, 0.5) * 2 / 1e-9)
ax2.set_xticklabels(list(map(str, numpy.arange(0.5, 3.5, 0.5))))
ax3.plot(rr, xi.flat[::-1])
# ax.set_xlim(10, 110)
ax.set_xlabel("Pore diameter (nm)")
ax.set_ylabel("$w(r_{\mathrm{G}})$")
ax2.set_xlabel("$\lambda_{\\mathrm{D}} / r_{\\mathrm{G}}$")
ax3.set_ylabel("$\\xi(r_{\\mathrm{G}})$")
ax3.set_ylim(0, 0.98)
# ax.set_xticks(numpy.arange(0.5, 3.5, 0.5) / 2 * lambda_D / 1e-9)
# plt.ylim(0, 0.04)
# plt.hist(repeat[::-1], xx[::-1])
# plt.ylim(0, 0.04)
 # plt.yticks([0, 0.02, 0.04])
plt.tight_layout()
plt.savefig(join(curdir, "../img/radius_dist_lamdbaD.svg"))
