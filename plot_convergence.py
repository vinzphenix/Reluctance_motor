import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

plt.rcParams['font.family'] = 'monospace'
fontsize = 14
fontsize2 = 13


if __name__ == "__main__":

    ############################################### - READ DATA - ###############################################

    filename = "data/convergence.txt"
    X = np.loadtxt(filename)
    nElem = X[:, 0]
    nNode = X[:, 1]
    L2_A = X[:, 2]
    H1_A = X[:, 3]
    lMIN = X[:, 4] * 1.e3
    lAVG = X[:, 5] * 1.e3
    lRMS = X[:, 6] * 1.e3
    lMAX = X[:, 7] * 1.e3

    L = lAVG

    filename = "data/convergenceNew.txt"
    X = np.loadtxt(filename)
    L2_B = X[:, 2]
    H1_B = X[:, 3]

    filename = "data/convergenceBis.txt"
    X = np.loadtxt(filename)
    L2_C = X[:, 2]
    H1_C = X[:, 3]

    ############################################### - GET ORDER - ###############################################

    model = LinearRegression().fit(
        np.log(L[:-1]).reshape((-1, 1)), np.log(L2_C))
    r_sq = model.score(np.log(L[:-1]).reshape((-1, 1)), np.log(L2_C))
    intercept = model.intercept_
    orderM, = model.coef_
    print("L2 order in n = {:.3f} \t R2 = {:.1f} \t intercept = {:.2f}".format(
        orderM, 100*r_sq, intercept))

    model = LinearRegression().fit(
        np.log(L[:-1]).reshape((-1, 1)), np.log(H1_C))
    r_sq = model.score(np.log(L[:-1]).reshape((-1, 1)), np.log(H1_C))
    intercept = model.intercept_
    orderM, = model.coef_
    print("H1 order in n = {:.3f} \t R2 = {:.1f} \t intercept = {:.2f}".format(
        orderM, 100*r_sq, intercept))

    ################################################ - FIGURE - #################################################

    fig, axs = plt.subplots(2, 1, figsize=(10, 6.5), constrained_layout='all')
    ax1 = axs[0]
    ax2 = axs[1]

    # FIRST PLOT
    ax1.loglog(L, L2_A, '-o', color='C0', label='$||\; a_i \;||$')
    ax1.loglog(L[:-1], L2_B, '-o', color='C1',
               label='$||\; a_{i} - a_{i+1} \;||$')
    ax1.loglog(L[:-1], L2_C, '-o', color='C2',
               label='$||\; a_{i} - a_{ref} \;||$')
    ax1.set_ylabel("norme $L^2$", fontsize=fontsize)

    # SECOND PLOT
    #ax2.loglog(nNode[1:], L2_B, '-o', color='C0', label='norme L2');
    ax2.loglog(L, H1_A, '-o', color='C0', label='$||\; a_i \;||$')
    ax2.loglog(L[:-1], H1_B, '-o', color='C1',
               label='$||\; a_{i} - a_{i+1} \;||$')
    ax2.loglog(L[:-1], H1_C, '-o', color='C2',
               label='$||\; a_{i} - a_{ref} \;||$')
    ax2.set_ylabel("norme $H^1$", fontsize=fontsize)

    ax2.set_xlabel("longueur caracteristique [mm]", fontsize=fontsize)

    for ax in axs:
        ax.grid(ls=':', which='both')
        ax.legend(fontsize=fontsize2)

    plt.savefig("./figures/convergence.svg", format='svg')
    plt.show()
