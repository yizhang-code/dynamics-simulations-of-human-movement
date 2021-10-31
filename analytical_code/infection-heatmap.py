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
sigma_r, sigma_theta = 2, 45 * np.pi / 180
R0 = 2
gamma = 1 / (14 * 24 * 3600.0)
beta_0 = gamma * R0
Nc = 10
p_daily = 15 / (24 * 60.0)
ρ_daily = Nc / (np.pi * 2**2) * p_daily
beta_max = beta_0 / (ρ_daily * sigma_r**2 * sigma_theta**2)

class t_sta(object):

    def __init__(self):
        self.count = np.zeros((max_i, max_i))
        self.total = 0
        self.dtheta = dtheta

    def add(self, theta1, theta2, val):
        theta1 += 180
        theta2 += 180
        i1, i2 = int(theta1 / self.dtheta), int(theta2 / self.dtheta)
        if i1 >= max_i or i2 >= max_i:
            return
        self.count[i1][i2] += val
        self.total += val

    def P(self, i1, i2):
        return self.count[i1][i2] * 120 * 120 / self.total


sta0 = t_sta()
files = [_ for _ in os.listdir(cur_dir + '/' + tar_dir) if '.dat' in _]
for file in files:
    print(file)
    with open('{}/{}/{}'.format(cur_dir, tar_dir, file), 'r') as rf:
        for i, line in enumerate(rf):
            r, agl1, agl2 = line.strip().split('\t')
            r, agl1, agl2 = float(r), float(agl1), float(agl2)
            if r <= 1.5 and r >= 0.25:
                beta = beta_max * np.exp(-(r ** 2 / (2 * sigma_r ** 2)) - (
                    (agl1 * np.pi / 180) ** 2 + (agl2 * np.pi / 180) ** 2) / (2 * sigma_theta ** 2))

                sta0.add(agl1, agl2, beta)
                if step == 2:
                    sta0.add(agl2, agl1, beta)

angle_arr = np.zeros((120, 120))
with open('{}/{}_infection_heatmap.dat'.format(cur_dir, tar_dir), 'w') as wf:
    for i in range(120):
        temp = []
        for j in range(120):
            temp.append(sta0.P(i, j).item())

        wf.write('\t'.join(str(ele) for ele in temp))
        wf.write('\n')
