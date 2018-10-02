import numpy
import matplotlib.pyplot as plt

file_name = "../data/pore-dist.csv"
data = numpy.genfromtxt(file_name, delimiter=",")

fig = plt.style.use("science")
plt.figure(figsize=(2.5, 2.5), facecolor="w")
repeat = [d for d, n in data for _ in range(int(n))]
plt.hist(repeat, numpy.arange(5,66,5), density=True)
plt.ylim(0, 0.04)
plt.yticks([0, 0.02, 0.04])
plt.xlabel("Pore Radius (nm)")
plt.ylabel("Frequencies")
plt.savefig("../img/radius_dist.svg")
