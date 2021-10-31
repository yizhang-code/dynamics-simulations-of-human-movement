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

tar_dir = sys.argv[1]
if tar_dir == 'oriAB':
    step = 1
else:
    step = 2
# tar_dir = 'ori'
dtheta, max_i = 3, 121
class t_sta0(object):

    def __init__(self):
        self.count = np.zeros((max_i, 1))
        self.total = 0
        self.dtheta = dtheta

    def add(self, theta):
        i = int(theta / self.dtheta)
        if i >= max_i:
            return
        self.count[i] += 1
        self.total += 1

    def P(self, i):
        return self.count[i] / self.total


class t_sta(object):

    def __init__(self):
        self.count = np.zeros((max_i, max_i))
        self.total = 0
        self.dtheta = dtheta

    def add(self, theta1, theta2):
        i1, i2 = int(theta1 / self.dtheta), int(theta2 / self.dtheta)
        if i1 >= max_i or i2 >= max_i:
            return
        self.count[i1][i2] += 1
        self.total += 1

    def P(self, i1, i2):
        return self.count[i1][i2] / self.total


sta0, sta_theta1, sta_theta2 = t_sta(), t_sta0(), t_sta0()
files = [_ for _ in os.listdir(cur_dir + '/' + tar_dir) if '.dat' in _]
for file in files:
    print(file)
    with open('{}/{}/{}'.format(cur_dir, tar_dir, file), 'r') as rf:
        for i, line in enumerate(rf):
            r, agl1, agl2 = line.strip().split('\t')
            r, agl1, agl2 = float(r), float(agl1) + 180, float(agl2) + 180
            if r <= 1.5 and r >= 0.25:
                sta0.add(agl1, agl2)
                sta_theta1.add(agl1)
                sta_theta2.add(agl2)
                if step == 2:
                    sta0.add(agl2, agl1)
                    sta_theta1.add(agl2)
                    sta_theta2.add(agl1)
angle_arr = np.zeros((120, 120))
with open('{}/{}_heatmap.dat'.format(cur_dir, tar_dir), 'w') as wf:
    for i in range(120):
        temp = []
        for j in range(120):
            if sta_theta1.P(i) != 0 and sta_theta2.P(j) != 0:
                # wf.write(sta0.P(i, j) / sta_theta1.P(i) / sta_theta2.P(j))
                temp.append(
                    (sta0.P(i, j) / sta_theta1.P(i) / sta_theta2.P(j)).item())
                # print(sta0.P(i, j) / sta_theta1.P(i) / sta_theta2.P(j))
            else:
                temp.append(0)
                # wf.write(0)
        wf.write('\t'.join(str(ele) for ele in temp))
        wf.write('\n')
