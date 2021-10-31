import numpy as np
import sys
import os

cur_dir = os.getcwd()

tar_dir = sys.argv[1]

N = 1000
dr = 0.05
r_cut = 0.01
max_r = 8


class Ave(object):
    def __init__(self):
        self.data = np.zeros((N, 2))

    def add(self, r, gr):
        indx_r = int((r - r_cut) / dr)
        if indx_r < 0 or indx_r > N:
            return
        self.data[indx_r][0] += gr
        self.data[indx_r][1] += 1

    def val(self, indx):
        if self.data[indx][1] == 0:
            return 0
        return self.data[indx][0] / self.data[indx][1]


ave_gr = Ave()

files = [_ for _ in os.listdir('{}/{}'.format(cur_dir, tar_dir)) if 'dat' in _]
for file in files:
    with open('{}/{}/{}'.format(cur_dir, tar_dir, file), 'r') as rf:
        for line in rf:
            line = line.split('\t')
            ave_gr.add(float(line[0]), float(line[1]))

with open('{}/{}_gr.csv'.format(cur_dir, tar_dir), 'w') as wf:
    for i in range(N):
        r = i * dr + r_cut
        ave_val = ave_gr.val(i)
        if ave_val != 0 and r < max_r:
            # print(i, i * dr + r_cut)
            wf.write('{},{}'.format(r, ave_val))
            wf.write('\n')
