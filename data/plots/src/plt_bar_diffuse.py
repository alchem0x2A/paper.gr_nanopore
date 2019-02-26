import numpy
import matplotlib.pyplot as plt


salts = ["KCl", "NaCl", "LiCl", "CaCl2", "MgSO4", "K2SO4", "KFeCN"]
data = numpy.genfromtxt("../data/diffuse-pcte-salt.csv",
                        delimiter=",",
                        skip_header=2)

plt.figure(figsize=(3, 3 * 0.8))
plt.style.use("science")

# for i in range(len(salts)):
    # name = salts[i]
    # Measured
data = data / 1e-6
print(data.shape)
x = numpy.arange(len(salts))
w = 1.0 /6
#bare pcte measure
plt.bar(x - w * 3 / 2 , data[:, 6], width=w,
        yerr=data[:, 7], color="blue",
        label="raw")

# Simulate
plt.bar(x - w / 2, data[:, 3], width=w, color="red",
        yerr=data[:, 5 : 3 : -1].T,
        label="simu")

# Bare
plt.bar(x + 1 / 2 * w, data[:, 1], width=w,
        yerr=data[:, 2], color="green",
        label="gr")



    # 
# plt.plot(data[:, 1] / 3600, data[:, 2] / 1e-3)
# plt.xlabel("t (h)")
# plt.ylabel("Conductivity")
plt.xticks(range(0, 7))
plt.xlim(-0.5, 6.5)

plt.legend()
plt.savefig("../img/diffusion-salts.svg")
