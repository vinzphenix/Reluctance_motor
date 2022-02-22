import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['font.family'] = 'monospace'
fontsize = 14
fontsize2 = 13


############################################### - READ DATA - ###############################################

filename = "data/torque.txt"
names = np.array(['coil', 'nonLinear', 'theta', 'torque'])
dfG = pd.read_csv(filename, delimiter=" ", names=names, index_col=False)
dfG['torque'] *= -1.

alpha = 0.5
phi0 = 18.33
dphi = 30.
dth = 0.5
limits = phi0 + np.arange(0., 181., 30.)
limits = np.roll(limits, 1)
limits[0] = 0.

fig, ax =  plt.subplots(1, 1, figsize=(10, 4.5), constrained_layout=True, sharex='all', sharey='all')
labels = ['bobines A', 'bobines B', 'bobines C']


############################################## - FIGURE Linear - ############################################

#ax = axs[0]
df = dfG[dfG['nonLinear'] == 0]
matrix = pd.pivot_table(df, values="torque", index=['theta'], columns=['coil'], aggfunc=np.mean)
theta = matrix.index
coils = matrix.columns
matrix = np.asmatrix(matrix)


idxs = [np.argwhere(
        (matrix[:, i] > matrix[:, (i+1)%3]) &(matrix[:, i] > matrix[:, (i+2)%3])
        )[:,0] for i in range(3)]

segments = {0: [], 1: [], 2: []}
for i in range(3):
    st = idxs[i][0]
    for j in range(len(idxs[i]) - 1):
        if idxs[i][j] + 1 != idxs[i][j+1]:
            segments[i].append((st, idxs[i][j]+1))
            st = idxs[i][j+1]
    segments[i].append((st, idxs[i][len(idxs[i]) - 1]+1))

sumTorque = 0.

for i in range(len(coils)):
    ax.plot(theta, matrix[:, i], color='C'+str(i), label=labels[i], alpha=0.5)
    for st, end in segments[i]:
        ax.plot(theta[st:end], matrix[st:end, i], color='C'+str(i), lw=2)
        sumTorque += np.sum(matrix[st:end, i])

mean = sumTorque / len(theta)
ax.axhline(mean, color='grey', ls='--', label='couple moyen')
ax.axvline(phi0, color='lightgrey', ls='--', label='18.33 °')

#ax.set_title("Relation Angle - Couple", fontsize=fontsize)
ax.set_ylabel("Couple [N m]", fontsize=fontsize)
ax.set_xlabel("Angle [deg]", fontsize=fontsize)


############################################ - FIGURE Non Linear - ##########################################
# if reactivated, code after needs to be modified
"""
ax = axs[1]
df = dfG[dfG['nonLinear'] == 1]

matrix = pd.pivot_table(df, values="torque", index=['theta'], columns=['coil'], aggfunc=np.mean)
theta = matrix.index
coils = matrix.columns
matrix = np.asmatrix(matrix)


idxs = [np.argwhere(
        (matrix[:, i] > matrix[:, (i+1)%3]) &(matrix[:, i] > matrix[:, (i+2)%3])
        )[:,0] for i in range(3)]

segments = {0: [], 1: [], 2: []}
for i in range(3):
    st = idxs[i][0]
    for j in range(len(idxs[i]) - 1):
        if idxs[i][j] + 1 != idxs[i][j+1]:
            segments[i].append((st, idxs[i][j]+1))
            st = idxs[i][j+1]
    segments[i].append((st, idxs[i][len(idxs[i]) - 1]+1))

sumTorque = 0.

for i in range(len(coils)):
    ax.plot(theta, matrix[:, i], color='C'+str(i), label=labels[i], alpha=0.5)
    for st, end in segments[i]:
        ax.plot(theta[st:end], matrix[st:end, i], color='C'+str(i), lw=2)
        sumTorque += np.sum(matrix[st:end, i])

mean = sumTorque / len(theta)
ax.axhline(mean, color='grey', ls='--')
ax.axvline(phi0, color='grey', ls='--')

ax.set_title("Couple - Non linéaire")
ax.set_ylabel("Couple [N m]")
ax.set_xlabel("Angle [deg]")
"""
################################################################################################

#for ax in axs :
ax.grid(ls=':')
ax.legend(loc='lower right', fontsize=fontsize2)

plt.savefig("./figures/torqueAngle.svg", format='svg')
plt.show()

