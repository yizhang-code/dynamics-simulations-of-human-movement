import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import os
import itertools
import seaborn as sns
from scipy.ndimage.filters import gaussian_filter
import matplotlib.patches as patches
import matplotlib as mpl

cur_dir = os.getcwd()

max_x, max_y = 73, 41

data = np.zeros((max_x, max_y))
with open('{}/r_theta.dat'.format(cur_dir), 'r') as rf:
    for i, line in enumerate(rf):
        line = line.strip().split('\t')
        for j, item in enumerate(line):
            data[j][i] = float(item)

fig, ax = plt.subplots(figsize=(9, 7), dpi=200)
array_smooth = gaussian_filter(data, sigma=1)
sns.heatmap(np.flip(array_smooth, axis=0),
            cmap="coolwarm", vmin=0, vmax=3.0, ax=ax)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)

levels = np.arange(1, 3, 1)
plt.contour(np.flip(array_smooth, axis=0),
            levels=levels, colors=['g', 'k'])

ax.set_ylabel(r'$\theta_1^{\degree}-\theta_2^{\degree}$', fontsize=45)
ax.set_xlabel('r(meters)', fontsize=45)
ax.set_xticks(np.linspace(0, 40, 5, dtype=int))
ax.set_xticklabels(np.linspace(0, 4, 5, dtype=int), fontsize=20)
ax.set_yticks(np.linspace(0, 72, 7, dtype=int))
ax.set_yticklabels(np.linspace(180, -180, 7, dtype=int), fontsize=20)
ax.set_xlim([0, 40])
ax.set_ylim([0, 72])
plt.gcf().subplots_adjust(bottom=0.15, left=0.18, top=0.95, right=0.99)
plt.show()
# plt.savefig('r_theta.png')
