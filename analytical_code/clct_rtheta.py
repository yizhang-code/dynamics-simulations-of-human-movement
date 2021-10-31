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

class t_sta(object):

    def __init__(self):
        self.heatmap = np.zeros((max_x, max_y))
        self.count = np.zeros((max_x, max_y))

    def add(self, i, j, val):
        if i >= max_y or j + 36 >= max_x:
            return
        self.heatmap[j + 36][i] += val
        self.heatmap[36 - j][i] += val
        self.count[j + 36][i] += 1
        self.count[36 - j][i] += 1

    def P(self, i, j):
        # for i in range(max_y):
        #     for j in range(max_x):
        count_ij = self.count[j][i]
        if count_ij != 0:
            return self.heatmap[j][i] / count_ij
        else:
            return 0
        # return self.heatmap


data = t_sta()

def load_file(file):

    with open('{}/r_theta_files/{}'.format(cur_dir, file), 'r') as rf:
        for i, line in enumerate(rf):
            line = line.strip().split(',')[:-1]
            for j, item in enumerate(line):
                if item != '':
                    data.add(i, j, float(item))


files = [_ for _ in os.listdir(
    '{}/r_theta_files'.format(cur_dir)) if 'dat' in _]


for file in files:
    load_file(file)

with open('{}/r_theta.dat'.format(cur_dir), 'w') as wf:
    for i in range(max_y):
        for j in range(max_x):
            wf.write('{}\t'.format(data.P(i, j)))
        wf.write('\n')

    # plt.figure(figsize=(9, 7))
    # fig, ax = plt.subplots(figsize=(9, 7))
    # array_smooth = gaussian_filter(data, sigma=1)
    # sns.heatmap(np.flip(array_smooth, axis=0),
    #             cmap="coolwarm", vmin=0, vmax=3.0, ax=ax)
    # cbar = ax.collections[0].colorbar
    # cbar.ax.tick_params(labelsize=20)

    # levels = np.arange(1, 3, 1)
    # plt.contour(np.flip(array_smooth, axis=0),
    #             levels=levels, colors=['g', 'k'])

    # ax.set_ylabel(r'$(\theta_1-\theta_2$)' + '\u00B0', fontsize=45)
    # ax.set_xlabel('r(meters)', fontsize=45)
    # ax.set_xticks(np.linspace(0, 40, 5, dtype=int))
    # ax.set_xticklabels(np.linspace(0, 4, 5, dtype=int), fontsize=20)
    # ax.set_yticks(np.linspace(0, 72, 7, dtype=int))
    # ax.set_yticklabels(np.linspace(180, -180, 7, dtype=int), fontsize=20)
    # ax.set_xlim([0, 40])
    # ax.set_ylim([0, 72])
    # plt.gcf().subplots_adjust(bottom=0.17, left=0.17, top=0.98, right=0.98)
    # # plt.show()
    # plt.savefig('{}_rtheta2.png'.format(file))
