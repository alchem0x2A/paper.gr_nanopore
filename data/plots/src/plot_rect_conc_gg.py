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

if __name__ == "__main__":
    delta0 = 2.2
    import sys
    try:
        name = sys.argv[1]
    except IndexError:
        name = "KCl"

    fig = plt.style.use("science")
    # plt.figure(figsize=(4, 4), facecolor="w")
    plt.figure(figsize=(2.8, 2.1), facecolor="w")
    # plt.axvline(x=0, ls="--", color="#00ffee")
    data, err = read(name)
    # data[-1, -1] *= 1.5
    # plt.errorbar(data[:, 0], data[:, 1], yerr=err / 2, fmt="o")
    conc = data[:, 0]
    eta = data[:, 1]
    xi = convert(delta=delta0, eta=eta)
    xi_err_lo = xi - convert(delta=delta0, eta=eta - err / 2)
    xi_err_hi = convert(delta=delta0, eta=eta + err / 2) - xi
    log_x = numpy.log10(conc)
    log_y = numpy.log10(xi)
    p = linregress(log_x, log_y)
    s = p.slope; b = p.intercept
    print(s, b, p.rvalue)
    plt.errorbar(conc, xi, yerr=(xi_err_lo, xi_err_hi), fmt="o")
    xx = numpy.logspace(-4.45, -1.5, 100)
    yy = 10 ** (numpy.log10(xx) * s + b)
    plt.plot(xx, yy, "--")
    # plt.xscale("log")
    plt.xlabel("KCl concentration (mol/L)")
    plt.ylabel("Rectification")
    plt.xscale("log")
    plt.yscale("log")
    plt.ylim(3e-2, 1)
    plt.savefig("../img/rect_{}_conc_gg.svg".format(name))
