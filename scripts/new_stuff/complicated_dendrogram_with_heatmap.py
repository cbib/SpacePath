import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import matplotlib.ticker as ticker

# Generate random features and distance matrix.
np.random.seed(10)  # for reproducible data
mydim = 2000
x = np.random.rand(mydim)
D = np.zeros([mydim, mydim])
for i in range(mydim):
    for j in range(mydim):
        D[i,j] = abs(x[i] - x[j])

condensedD = squareform(D)
Y = sch.linkage(condensedD, method='centroid')

## set linewidth
matplotlib.rcParams['lines.linewidth'] = 0.7

# Compute and plot left dendrogram.
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_axes([0.15, 0.1, 0.095, 0.6])

Z1 = sch.dendrogram(Y, orientation='left')
ax1.set_xticks([])
ax1.set_yticks([])

# Compute and plot top dendrogram.
ax2 = fig.add_axes([0.25, 0.71, 0.4, 0.2])
Z2 = sch.dendrogram(Y)
ax2.set_xticks([])
ax2.set_yticks([])

# Plot distance matrix.
axmatrix = fig.add_axes([0.25, 0.1, 0.4, 0.6])
idx1 = Z1['leaves']
idx2 = Z2['leaves']
D = D[idx1,:]
D = D[:,idx2]

mask = np.triu(np.ones_like(D))
# np.flip((1 - D) * mask, axis=0)
im = axmatrix.matshow(1 - D,
                      aspect='auto', origin='lower',
                      cmap=plt.cm.Blues)
labels_truncated = [f"uuuuu uuu uuu ' % 5 s [ 2343  ]u--{j}" for j in range(D.shape[0])]

ticks= [int(j) for j in range(D.shape[0])]
autolocator = int(D.shape[0]/40)
if autolocator <= 0:
    autolocator = 1
axmatrix.yaxis.set_major_locator(ticker.MultipleLocator(autolocator))
axmatrix.set_yticklabels(labels_truncated)  # put  #  , rotation=45
axmatrix.yaxis.set_label_position("right")
axmatrix.yaxis.tick_right()
axmatrix.set_xticks([])  # remove axis labels

# Plot colorbar.
#axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.6])
axcolor = fig.add_axes([0.05, 0.1, 0.015, 0.3 ])
plt.colorbar(im, cax=axcolor)
plt.show()

print("popop", bbox_inches='tight')

