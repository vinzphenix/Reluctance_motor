import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['font.family'] = 'monospace'
fontsize = 14
fontsize2 = 13


if __name__ == "__main__":

    Rin = np.loadtxt("./build/radiusIn.txt")
    Rout = np.loadtxt("./build/radiusOut.txt")

    n = np.arange(len(Rin))
    plt.hist(Rin, 500)

    n = np.arange(len(Rout))
    plt.hist(Rout, 500)

    plt.show()
