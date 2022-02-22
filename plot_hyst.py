import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt

plt.rcParams['font.family'] = 'monospace'
fontsize = 14
fontsize2 = 13


if __name__ == "__main__":

    ############################################### - READ DATA - ###############################################

    M = np.loadtxt("./data/mu.txt")

    b2 = M[:, 0]
    b = sqrt(b2)
    invMu = M[:, 1]
    dinvMu = M[:, 2]

    H0 = np.array([
        0, 10, 20, 30, 40, 50,
        60, 70, 80, 90, 100, 125, 150, 175, 200, 250,
        300, 400, 500, 600,  700, 800, 900, 1000, 1250, 1500, 2000, 2500, 5000,
        7500,  10000, 15000, 20000, 59000, 174000, 514000, 1520000, 4470000,
        13200000, 38900000, 115000000, 339000000, 1000000000])

    B0 = np.array([
        0.0,
        0.194880829963, 0.377143018857, 0.537767739762, 0.672888260835,
        0.783043000477, 0.871342430831, 0.941778611986, 0.998183303557, 1.04378111223,
        1.08110469369, 1.14963767549, 1.19607212343, 1.22964695907, 1.25515221835,
        1.29162498935, 1.31678879432, 1.35015120537, 1.37220092877, 1.38859114656,
        1.4017440574, 1.41287024565, 1.42264180514, 1.43146158921, 1.45082466146,
        1.46784549989, 1.49819370601, 1.52578650709, 1.64314027719, 1.73458485332,
        1.8039068939, 1.89568786291, 1.95213815187, 2.1390774927, 2.45827909293,
        3.32303272825, 5.85485500678, 13.2701832298, 35.2114648741, 99.8027446541,
        291.062951228, 854.036370229, 2515.3105707])

    end = 33  # len(H0)

    H = H0[:end]
    B = B0[:end]

    ################################################ - FIGURE - #################################################

    fig, axs = plt.subplots(1, 3, figsize=(10, 3.5), constrained_layout=True)

    ax = axs[0]
    ax.set_title("B( H )", fontsize=fontsize)
    ax.plot(H[:end], B[:end], 'o', label='data')
    ax.plot(b * invMu, b, '-', label='spline')
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    #ax.set_xlabel("H [A/m]")
    #ax.set_ylabel("B [T]")

    ax = axs[1]
    ax.set_title(r"$\mu^{-1}$ ( B )", fontsize=fontsize)
    #ax.plot(B[1:]*B[1:], H[1:]/B[1:], 'o', label='data')
    #ax.plot(b2, invMu, '-', label='spline')
    ax.plot(B[1:], H[1:]/B[1:], 'o', label='data')
    ax.plot(b, invMu, '-', label='spline')

    ax = axs[2]
    ax.set_title(r"$d (\mu^{-1}) / dB^2 $ ( B )", fontsize=fontsize)
    #ax.plot(B[1:]*B[1:], np.gradient(H0[1:]/B0[1:], B0[1:]*B0[1:])[0:end-1], 'o', label='data')
    #ax.plot(b2, dinvMu, '-', label='spline')
    ax.plot(B[1:], np.gradient(H0[1:]/B0[1:], B0[1:]*B0[1:])
            [0:end-1], 'o', label='data')
    ax.plot(b, dinvMu, '-', label='spline')

    for ax in axs:
        ax.grid(ls=':')
        ax.legend(fontsize=fontsize2)

    plt.savefig("./figures/hysteresis.svg", format='svg')
    plt.show()
