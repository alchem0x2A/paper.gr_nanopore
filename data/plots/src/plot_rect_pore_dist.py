from utils import *
import numpy
import matplotlib.pyplot as plt
import os
from os.path import dirname, join, exists
from scipy.constants import pi, hbar, e
import numpy
from scipy.interpolate import griddata, interp2d
import pickle

curdir = dirname(__file__)

"""
Plot rectification based on pore distribution
"""

func_rec_interp = pickle.load(open(join(curdir,
                                   "../data/FEM/concentration/1D",
                                   "rect_2d_intep.pickle"), "rb"))

data_pore = numpy.genfromtxt(join(curdir, "../data/exp/pore-dist.csv"),
                       delimiter=",")
r_exp = data_pore[:, 0]
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
