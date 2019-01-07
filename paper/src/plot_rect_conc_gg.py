import numpy
import matplotlib.pyplot as plt
from scipy.stats import linregress

file_templ = "../data/rect_{}_conc.csv"
def read(name):
    file_name = file_templ.format(name)
    err = None
    with open(file_name) as f:
        err = float(f.readline().strip())
    data = numpy.genfromtxt(file_name,
                            skip_header=1,
                            delimiter=",")
    return data, err

def convert(delta, eta):
    assert delta > 0
    xi_inv = (1 - eta) / (delta * eta + 1)
    return 1 - xi_inv

import sys
try:
    name = sys.argv[1]
except IndexError:
    name = "KCl"

fig = plt.style.use("science")
plt.figure(figsize=(4, 4), facecolor="w")
# plt.axvline(x=0, ls="--", color="#00ffee")
data, err = read(name)
# data[-1, -1] *= 1.5
# plt.errorbar(data[:, 0], data[:, 1], yerr=err / 2, fmt="o")
conc = data[:, 0]
eta = data[:, 1]
for delta in (0.2, 0.5, 1, 2, 5, 10, 20, 30):
    xx = conc
    yy = convert(delta, eta)
    log_x = numpy.log(xx)
    log_y = numpy.log(yy)
    p = linregress(log_x, log_y)
    s = p.slope
    plt.plot(conc, convert(delta, eta),
             "o",
             label="$\\delta={{{0}}}, p={{{1:.2f}}}$".format(delta, s))
# plt.plot(xx, yy, "--")
# plt.xscale("log")
plt.xlabel("KCl concentration (mol/L)")
plt.ylabel("Rectification")
plt.xscale("log")
plt.yscale("log")
plt.legend(loc=0)
plt.savefig("../img/rect_{}_conc_gg.svg".format(name))
