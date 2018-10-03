import numpy
import matplotlib.pyplot as plt

def get_time(string):
    string = string.decode("ascii")
    h, m, s = map(float, string.split(":"))
    return 3600 * h + 60 * m + 1 * s
    
data = numpy.genfromtxt("../data/diffusion-slope.csv",
                        converters={1: get_time},
                        delimiter=",")

plt.figure(figsize=(3, 3 * 0.8))
plt.style.use("science")
cond = numpy.where(data[:, 1] < 3600)
# plt.plot(data[:, 1][cond] / 3600, data[:, 2][cond] / 1e-3)
# plt.plot(data[:, 1][data[:, 1] > 7200] / 3600, data[:, 2][data[:, 1] > 7200] / 1e-3)
plt.plot(data[:, 1] / 3600, data[:, 2] / 1e-3)
plt.xlabel("t (h)")
plt.ylabel("Conductivity")
plt.savefig("../img/diffusion-slope.svg")
