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

fig = plt.style.use("science")
plt.figure(figsize=(2.8, 2.8), facecolor="w")
plt.axvline(x=0, ls="--", color="#00ffee")
plt.axhline(y=0, ls="--", color="#00ffee")

for name in ["KCl", "NaCl", "LiCl",
             "MgSO4", "K2SO4", "CaCl2",
             "K3FeCN"]:

    data, err = read(name)
    plt.errorbar(data[:, 0], data[:, 1], yerr=err / 2,
                 fmt="o-", alpha=0.5, label=name,
                 markersize=4)
    # plt.fill_between(data[:, 0], data[:, 1] - err / 2,
                     # data[:, 1] + err / 2,
                     # fmt="o",
                     # alpha=0.5, label=name)
plt.ylim(-0.1, 0.7)
plt.xlim(-1.25, 1.25)
plt.xlabel("$V_{\\mathrm{G}}$ (V)")
plt.ylabel("Rectification")
plt.legend()
plt.savefig("../img/rect_all.svg".format(name))
