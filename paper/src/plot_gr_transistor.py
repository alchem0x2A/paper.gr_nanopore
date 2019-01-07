import numpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.constants import epsilon_0
from scipy.constants import elementary_charge as e
from scipy.constants import hbar, pi

file_name = "../data/gr_sample.csv"
data = numpy.genfromtxt(file_name, delimiter=",")
# print(data)
VG = data[:, 0]
IDS = data[:, -1]
vv = numpy.linspace(VG.min(), VG.max(), 128)
ii = interp1d(VG, IDS, kind="linear")(vv)
p = numpy.polyfit(vv[(vv < 4) & (vv > 0)], ii[(vv < 4) & (vv > 0)], 2)
print(-p[1]/(2*p[0]))

v_min = -p[1]/(2*p[0])
print("Dirac point", v_min)

plt.figure(figsize=(3.2, 2.5))
plt.style.use("science")

plt.plot(VG, IDS / 1e-6, "o", markersize=3)
plt.xlabel("$V_{\\rm{G}}$ (V)")
plt.ylabel("$I_{\\rm{DS}}$ ($\\rm{\\mu}$)A")
plt.axvline(x=v_min, ls="--")
plt.xlim(-80, 80)
plt.tight_layout()
plt.savefig("../img/pg-transistor.svg")

epsilon_sio2 = 2.8
d = 300e-9
c = epsilon_0 * epsilon_sio2 / d
sigma = v_min * c
sigma_e = sigma / e / 1e4 / 1e12
print(sigma, sigma_e)
vf = 1.1e6
del_v = hbar * vf / e * numpy.sqrt(pi * abs(sigma) / e)
print(del_v)
