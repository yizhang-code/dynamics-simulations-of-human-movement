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

tar_file = sys.argv[1]

max_i = 120
data = np.zeros((max_i, max_i))
with open('{}/{}_infection_heatmap.dat'.format(cur_dir, tar_file), 'r') as rf:
    for i, line in enumerate(rf):
        line = line.strip().split('\t')
        for j, item in enumerate(line):
            data[i][j] = float(item)

plt.figure(figsize=(9, 7))
array_smooth = gaussian_filter(data, sigma=1)
ax = sns.heatmap(np.flip(array_smooth, axis=0),
                 cmap="coolwarm", vmin=0.0, vmax=10)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)

levels = np.arange(1, 2, 1)
print(levels)
plt.contour(np.flip(array_smooth, axis=0), levels=levels, colors=['y', 'k'])

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
theta1 = '\u03F4' + '1'
theta2 = '\u03F4' + '2'
s1 = theta1.translate(SUB)
s2 = theta2.translate(SUB)
if tar_file == 'ori':
    ylabel = s2 + '\u00B0' + '/Children'
    xlabel = s1 + '\u00B0' + '/Children'
elif tar_file == 'ori_ppl':
    ylabel = s2 + '\u00B0' + '/individuals'
    xlabel = s1 + '\u00B0' + '/individuals'
else:
    ylabel = s2 + '\u00B0' + '/Teachers'
    xlabel = s1 + '\u00B0' + '/Children'

plt.ylabel(ylabel, fontsize=45)
plt.xlabel(xlabel, fontsize=45)
plt.xticks(np.arange(0, 121, 20), np.arange(-180,
                                            181, 60), fontsize=20, rotation=0)
plt.yticks(np.arange(0, 121, 20), np.arange(
    180, -181, -60), fontsize=20, rotation=0)
plt.gcf().subplots_adjust(bottom=0.17)
plt.gcf().subplots_adjust(left=0.17)
plt.gcf().subplots_adjust(top=0.98)
plt.gcf().subplots_adjust(right=0.98)
plt.show()
