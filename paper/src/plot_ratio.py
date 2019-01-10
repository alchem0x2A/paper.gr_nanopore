import numpy
import matplotlib.pyplot as plt

file_templ = "../data/{}.csv"
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

def eta_to_xi(eta, delta=2.20):
    return (delta + 1) * eta / (delta * eta + 1)

def main():
    import sys
    try:
        name = sys.argv[1]
    except IndexError:
        name = "KCl"

    fig = plt.style.use("science")
    plt.figure(figsize=(2.8, 2.1), facecolor="w")
    plt.axvline(x=0, ls="--", color="#00ffee")
    plt.axhline(y=0, ls="--", color="#00ff00")
    data, err = read(name)
    xi = eta_to_xi(data[:, 1])
    err_lo = xi - eta_to_xi(data[:, 1] - err / 2)
    err_hi = eta_to_xi(data[:, 1] + err / 2) - xi
    plt.errorbar(data[:, 0], eta_to_xi(data[:, 1]), yerr=(err_lo, err_hi), fmt="o")
    plt.xlim(-1.3, 1.3)
    plt.ylim(-0.15, 0.6)
    plt.xlabel("$V_{\mathrm{G}}$ (V)")
    plt.ylabel("Rectification")
    plt.savefig("../img/rect_{}.svg".format(name))

if __name__ == "__main__":
    main()
