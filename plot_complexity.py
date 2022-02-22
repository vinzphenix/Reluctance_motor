import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression
from matplotlib import cm

cmap = cm.get_cmap('viridis_r')
fontsize = 14
fontsize2 = 13


if __name__ == "__main__":

    ############################################### - READ DATA - ###############################################

    filename = './data/complexity.txt'
    df = pd.read_csv(filename, delimiter=" ", header=0, index_col=False)
    df.rename(columns={'matrixVectorProduct': 'A * d_k'}, inplace=True)

    df.drop(columns=["PCG"], inplace=True)
    #df.drop(columns=['matrixVectorProduct', 'solve', 'ILU', 'sort'], inplace=True)

    s = df.columns.drop(["nNode", "nIter"])
    df.loc[:, s] = df.loc[:, s].mul(1./df['nIter'], 0)

    df.insert(len(df.columns)-1, 'other',
              df['full'] - df[s.drop("full")].sum(axis=1), True)
    s = df.columns.drop(["nNode", "nIter"])

    df = df.groupby(['nNode'], as_index=False).mean()

    #matrix = 1.e3 * pd.pivot_table(df, values="time", index=['theta'], columns=['nNode'], aggfunc=np.mean)
    matrix = np.zeros((len(df.nNode), 1+len(s)))
    matrix[:, 1:] = np.array(df[s])
    #matrix = np.cumsum(matrix, axis=1)

    ############################################### - GET ORDER - ###############################################

    X = np.log(np.array(df['nNode'])).reshape((-1, 1))
    Y = np.log(np.array(df['full']))
    model = LinearRegression().fit(X, Y)
    order, = model.coef_
    intercept = model.intercept_

    ################################################ - FIGURE - #################################################

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    # FIRST PLOT
    ax = axs[0]
    ax.loglog(df['nNode'], 1e3 * df['full'], '-o', label='full time')
    xx = np.linspace(np.array(df['nNode'])[0], np.array(df['nNode'])[-1])
    ax.loglog(xx, 1e3*np.exp(intercept) * (xx**order), ls='--',
              label='order {:.2f}'.format(order), alpha=0.75)

    ax.grid(ls=':')
    ax.set_xlabel("number of nodes", size=fontsize)
    ax.set_ylabel("time [ms]", size=fontsize)
    ax.legend(fontsize=fontsize2)

    # SECOND PLOT
    ax = axs[1]
    colors = [cmap(0. + i/7 * 1.) for i in range(6)] + \
        ['C3', 'C1', 'C0', 'lightgrey']
    width = 0.25

    ind = np.arange(5)
    for i in range(len(s) - 1):
        y1 = matrix[:, i+1] / matrix[:, -1]
        #y2 = matrix[:, i] / matrix[:, -1]
        #ax.fill_between(df['nNode'], y1, y2, label=df.columns[i+1])
        bottom = 1e2 * np.sum(matrix[:, :i+1], axis=1) / matrix[:, -1]
        ax.bar(ind, 1e2 * y1, width, bottom=bottom,
               label=s[i], color=colors[i])

    ax.set_xticks(ind)
    ax.set_xticklabels(np.array(df.nNode))

    ax.grid(ls=':')
    ax.set_xlabel("number of nodes", size=fontsize)
    ax.set_ylabel("répartition [%]", size=fontsize)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1],
              prop={'family': 'monospace', 'size': fontsize2}, ncol=2)

    plt.savefig("./figures/complexity.svg", format='svg')
    plt.show()
