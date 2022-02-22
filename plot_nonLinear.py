import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter


plt.rcParams['font.family'] = 'monospace'
fontsize = 14
fontsize2 = 13

if __name__ == "__main__":
    ############################################### - READ DATA - ###############################################

    filename = "data/nonLinear.txt"
    names = np.array(['nNode', 'flag', 'bobine', 'js', 'theta', 'torque'])
    dfG = pd.read_csv(filename, delimiter=" ",  names=names, index_col=False)
    dfG['torque'] *= -1.
    dfG['theta'] = -np.degrees(dfG['theta'])

    js0 = 8.8464 * 1.e5
    dfG = dfG[np.abs(dfG['js'] - js0) > 1.]
    dfG = dfG[np.abs(dfG['js'] - 81.*js0) > 1.]

    df1 = dfG[dfG['flag'] == 0]
    df2 = dfG[dfG['flag'] == 1]

    fig = plt.figure(figsize=(12, 9), constrained_layout='all')
    gs = GridSpec(5, 2, figure=fig)

    ############################################ - FIGURE Non Linear - ##########################################

    fontdict = {'fontsize': 13, 'fontweight': 'normal',
                'fontfamily': 'monospace'}

    matrix1 = pd.pivot_table(df1, values="torque", index=[
                             'theta'], columns=['js'], aggfunc=np.mean)
    matrix2 = pd.pivot_table(df2, values="torque", index=[
                             'theta'], columns=['js'], aggfunc=np.mean)
    theta = np.array(matrix1.index)
    js = np.array(matrix1.columns)

    matrix1 = np.array(matrix1)
    matrix2 = np.array(matrix2)

    for i in range(len(js)):
        if i == 0:
            ax = fig.add_subplot(gs[i, 0])
        else:
            ax = fig.add_subplot(gs[i, 0], sharex=ax)

        ax.text(0.85, 0.83, "k = {:2.0f}".format(js[i]/js0), fontsize=fontsize,
                bbox=dict(facecolor='forestgreen',
                          edgecolor='none', alpha=0.5),
                transform=ax.transAxes)
        ax.plot(theta, matrix1[:, i], color='C0', ls='--', label="linear")
        ax.plot(theta, matrix2[:, i], color='C1', ls='-', label="nonlinear")

        ax.set_ylabel(" ", fontsize=fontsize)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    ax.set_xlabel("Angle [degré]", fontsize=fontsize)

    ################################################# - ORDER - #################################################

    X = np.log(js).reshape((-1, 1))
    Y = np.log(np.mean(matrix2, axis=0))
    model = LinearRegression().fit(X, Y)
    order, = model.coef_
    intercept = model.intercept_
    print(order)

    ################################################################################################

    ax = fig.add_subplot(gs[:, 1])
    y1 = np.mean(matrix1, axis=0)
    y2 = np.mean(matrix2, axis=0)
    ax.loglog(js, y1, '-o', ls='--', label='solution linéaire')
    ax.loglog(js, y2, '-o', ls='-', label='solution non linéaire')

    ax.set_xlabel("densité de courant [A $s^{-2}$]", fontsize=fontsize)
    ax.set_ylabel("Couple [N m]", fontsize=fontsize)
    ax.set_title("Couple moyen", fontsize=fontsize+1)

    for ax in fig.get_axes():
        ax.grid(ls=':')

    ax.legend(loc='lower right', fontsize=fontsize2)
    fig.text(0.01, 0.5, 'Couple [N m]', va='center',
             rotation='vertical', fontsize=fontsize)

    (fig.get_axes()[0]).set_title(
        "Relation Angle - Couple", fontsize=fontsize+1)
    # Impact de la non-linéarité sur le couple

    plt.savefig("./figures/nonLinear.svg", format='svg')
    plt.show()
