import numpy
from scipy.constants import pi, hbar, e
vf = 1.1e6

# Use C/m^2
def delta_phi_gr(sigma):
    fac = hbar * vf / e * numpy.sqrt(pi * numpy.abs(sigma) / e)
    return fac * numpy.sign(sigma)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    data = numpy.genfromtxt("../result/radius/10/1D/sigma_0.001.txt", comments="%")
    phi_surf = data[:, 0]
    sigma_gr = -data[:, 1]
    d_phi  = delta_phi_gr(sigma_gr)
    V_b = phi_surf + d_phi
    fig = plt.figure(figsize=(2.8, 2.1))
    plt.style.use("science")
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(V_b, sigma_gr, "s-", markersize=4)
    ax2.plot(V_b, phi_surf, "o-", color="r", markersize=4)
    
    ax2.set_ylabel(r"$\phi_{\rm{0}}$ (V)")
    ax1.set_ylabel(r"$\sigma_{\rm{G}}$ (C/m$^{2}$)")

    ax1.set_xlabel("$V_{\\mathrm{G}}$ (V)")
    plt.xlim(0, 1.25)
    ax1.set_ylim(0, 0.1)
    ax2.set_ylim(0, 0.3)
    # plt.axhline(y=1.25, color="red", alpha=0.6, ls="--")

    plt.tight_layout()
    
    plt.savefig("two_panels.svg")
