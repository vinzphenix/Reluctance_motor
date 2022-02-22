import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.lines import Line2D
from scipy.sparse import coo_matrix
from scipy.interpolate import splev, splrep
from matplotlib.patches import Rectangle

fontsize = 14
fontsize2 = 13


############################################### - READ DATA - ###############################################

fileFULL = './data/renumber2FULL.txt'
fileONCE = './data/renumber2ONCE.txt'
fileZERO = './data/renumber2ZERO.txt'
pd.options.display.max_rows = 100

files = [fileFULL, fileONCE, fileZERO]
names = np.array(['nNode', 'theta', 'time'])

#df1 = pd.read_csv(fileFULL, delimiter=" ", names=names, index_col=False)
#df1 = df1.groupby(['nNode', 'theta'], as_index=False).mean()

dfs = [pd.read_csv(filename, delimiter=" ", names=names, index_col=False) for filename in files]

for df in dfs:
    df.drop(df[df.theta < 1e-5].index, inplace=True)
    df.loc[df.theta > 2 * np.pi - 1e-3, 'theta'] = 0.


matrices = [1.e3 * pd.pivot_table(df, values="time", index=['theta'], columns=['nNode'], aggfunc=np.mean) for df in dfs]
nNode = matrices[0].columns


############################################## - FIGURE TIME - ##############################################

fig, ax =  plt.subplots(1, 1, figsize=(6, 8), constrained_layout=True, sharex='all', sharey='all')

for i, size in enumerate(nNode):
    line = ax.plot(np.degrees(matrices[0].index), matrices[0][size], '-o', color='C'+str(i), ls= '-', label="{:>6d} noeuds".format(size))
    line = ax.plot(np.degrees(matrices[1].index), matrices[1][size], '-o', color='C'+str(i), ls='--')#, label="{:d} once".format(size))
    line = ax.plot(np.degrees(matrices[2].index), matrices[2][size], '-o', color='C'+str(i), ls=':')#, label="{:d} none".format(size))
    ax.set_yscale("log")


handles, labels = ax.get_legend_handles_labels()
first_legend = ax.legend(handles[::-1], labels[::-1], loc='center right', prop={'family': 'monospace', 'size': fontsize2})

legend_elements = [Line2D([0], [0], color='k', ls=':' , lw=2, label='None'),
                   Line2D([0], [0], color='k', ls='--', lw=2, label='Once'),
                   Line2D([0], [0], color='k', ls='-' , lw=2, label='Full')]
ax2 = ax.add_artist(first_legend)
plt.legend(handles=legend_elements, loc='center left', fontsize=fontsize2)

ax.grid(ls=':')
ax.set_xlabel(r"$\theta$ [°]", fontsize=fontsize)
ax.set_ylabel("time [ms]", fontsize=fontsize)

plt.savefig("./figures/numberingTime1.svg", format='svg')
plt.show()


############################################# - FIGURE SPARSE - #############################################
#SLOW -> in comments
"""
#A = np.loadtxt('./data/cooMatrixFULL.txt')
A = np.loadtxt('./data/cooMatrixONCE.txt') # most interresting one
#A = np.loadtxt('./data/cooMatrixZERO.txt')
print("loaded")
fig, axs = plt.subplots(4, 3, figsize=(6, 8), constrained_layout=True, sharex='all', sharey='all')


s = 0
for i, ax in enumerate(axs.flatten()):
    
    n = int(A[s, 0])
    nnz = int(A[s, 1])
    th = A[s, 2]
    
    row = A[s+1 : s+1+nnz, 0]
    col = A[s+1 : s+1+nnz, 1]
    data = A[s+1 : s+1+nnz, 2]
    
    sparse = coo_matrix((data, (row, col)), shape=(n, n))
    
    title = "{:.0f}°".format(np.degrees(th))
    if np.degrees(th) < 10 :
        pos = 0.92
    elif np.degrees(th) < 100:
        pos = 0.89
    else:
        pos = 0.85
    text = ax.text(pos,.89, title, size=12, horizontalalignment='center', transform=ax.transAxes)
    text.set_bbox(dict(facecolor='wheat', alpha=0.5, edgecolor='none'))

        
    ax.spy(sparse, markersize=1, rasterized=True)
    
    s += nnz + 1
    
    ax.set_yticklabels([])
    ax.set_xticklabels([])

plt.savefig("./figures/numberingSpy.svg", format='svg', dpi=300)
plt.show()"""

