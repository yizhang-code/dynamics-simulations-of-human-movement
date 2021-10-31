import os
import sys
# import random
from math import sqrt
import numpy as np
# import matplotlib.pyplot as plt

cur_dir = os.getcwd()
dt = float(sys.argv[1])
max_t = int(sys.argv[1])


class LoadAllxy(object):
    """docstring for LoadAllxy"""

    def __init__(self, particle_indx: int):
        self.loc = []
        self.particle_indx = particle_indx
        self.particle_num = 0

    @staticmethod
    def exist(x, y):
        if x != -1 and y != -1:
            return True
        else:
            return False

    @staticmethod
    def find_center(lx, ly, rx, ry):
        lx, ly, rx, ry = float(lx), float(ly), float(rx), float(ry)
        if LoadAllxy.exist(lx, ly) and LoadAllxy.exist(rx, ry):
            cx, cy = (lx + rx) / 2, (ly + ry) / 2
            dx, dy = rx - lx, ry - ly
            ds = sqrt(dx**2 + dy**2)
            if ds != 0:
                flag = 1
                Tx, Ty = -dy / ds, dx / ds
                abs_ori = np.arctan2(Ty, Tx)
            else:
                flag = -1
                abs_ori = None
        else:
            flag = -1
            cx, cy = -1, -1
            abs_ori = None
        return (flag, cx, cy, abs_ori)

    def extract(self, sec_loc):
        lx, ly, rx, ry = sec_loc[0], sec_loc[1], sec_loc[2], sec_loc[3]
        c_xy = LoadAllxy.find_center(lx, ly, rx, ry)
        # print(lx, ly, rx, ry, c_xy)
        self.loc.append(c_xy)

    def load_all_xy(self):
        with open('{}/{}'.format(cur_dir, 'all_xy.csv'), 'r') as rf:
            for _, line in enumerate(rf):
                line = line.split(',')[1:]
                self.extract(line[self.particle_indx *
                                  4:self.particle_indx * 4 + 4])
            self.particle_num = int(len(line) / 4)
        return self.loc


class ClctMSD(object):
    """docstring for ClctMSD"""

    def __init__(self, loc, data: dict, particle_num: int):
        self.loc = loc
        self.max_t = max_t
        self.data = data
        self.chunk_loc = {}
        self.chunk_ori = {}
        self.particle_num = particle_num

    def chunk_raw(self):
        flag = 0
        self.chunk_loc[flag], self.chunk_ori[flag] = [], []
        for i, item in enumerate(self.loc):
            if item[0] != -1:
                self.chunk_loc[flag].append((item[1], item[2]))
                self.chunk_ori[flag].append(item[3])
            if len(self.chunk_loc[flag]) != 0 and item[0] == -1:
                flag += 1
                self.chunk_loc[flag] = []
                self.chunk_ori[flag] = []
            else:
                pass

    @staticmethod
    def clct_dis(x1, y1, x2, y2):
        return (x1 - x2)**2 + (y1 - y2)**2

    @staticmethod
    def map2range(ori_list):
        if len(ori_list) == 0:
            return ori_list
        ranged_ori_list = [ori_list[0]]
        for i in range(1, len(ori_list)):
            increment = ori_list[i] - ori_list[i - 1]
            if increment > np.pi:
                increment -= 2 * np.pi
            if increment < -np.pi:
                increment += 2 * np.pi
            ranged_ori_list.append(ranged_ori_list[i - 1] + increment)
        return ranged_ori_list

    def output(self, ori_list, flag):
        if not os.path.exists('{}/abs-angles'.format(cur_dir)):
            os.mkdir('{}/abs-angles'.format(cur_dir))
        if len(ori_list) <= 600:
            return True
        else:
            with open('{}/abs-angles/{}_{}_ori.csv'.format(cur_dir, self.particle_num, flag), 'w') as wf:
                # with open('{}/abs-angles/{}_ori.csv'.format(cur_dir, 'abs-agl-test',  flag), 'w') as wf:
                for ori in ori_list:
                    wf.write('{}\n'.format(ori))

    def clctmsd(self, t):
        sum_dx, sum_dy, sum_da = 0, 0, 0
        sum_dr2, sum_da2, N = 0, 0, 0
        sum_dxda, sum_dyda, sum_dxdy = 0, 0, 0
        for flag, loc_list in self.chunk_loc.items():
            ori_list = self.map2range(self.chunk_ori[flag])
            for i in range(1, len(ori_list)):
                if ori_list[i] - ori_list[i - 1] > np.pi:
                    print(i, ori_list[i], ori_list[i - 1])
            self.output(ori_list, flag)
            if len(loc_list) > t:
                for i in range(len(loc_list) - t):
                    dx = loc_list[i + t][0] - loc_list[i][0]
                    dy = loc_list[i + t][1] - loc_list[i][1]
                    da = ori_list[i + t] - ori_list[i]

                    sum_dr2 += dx * dx + dy * dy
                    sum_da2 += da * da
                    sum_dx += dx
                    sum_dy += dy
                    sum_da += da
                    sum_dxdy += dx * dy
                    sum_dxda += dx * da
                    sum_dyda += dy * da
                    N += 1
        if N == 0:
            return (0, 0, 0, 0, 0, 0, 0, 0, 0)
        else:
            return (sum_dr2 / N, sum_da2 / N, sum_dx / N, sum_dy / N, sum_da / N, sum_dxdy / N, sum_dxda / N, sum_dyda / N, 1)

    def find_msd_vs_t(self):
        self.chunk_raw()
        for t in range(self.max_t):
            print(self.particle_num, t)
            dr2, da2, dx, dy, da, dxdy, dxda, dyda, count = self.clctmsd(t)
            if t not in self.data:
                self.data[t] = [0, 0, 0, 0, 0, 0, 0, 0, 0]
            self.data[t][0] += dr2
            self.data[t][1] += da2
            self.data[t][2] += dx
            self.data[t][3] += dy
            self.data[t][4] += da
            self.data[t][5] += dxdy
            self.data[t][6] += dxda
            self.data[t][7] += dyda
            self.data[t][8] += count
        return self.data


if __name__ == '__main__':
    data = {}
    particle_0 = LoadAllxy(0)
    particle_0_loc = particle_0.load_all_xy()
    particle_num = particle_0.particle_num
    data = ClctMSD(particle_0_loc, data, 0).find_msd_vs_t()

    for num in range(1, particle_num):
        print(num)
        data = ClctMSD(LoadAllxy(num).load_all_xy(), data, num).find_msd_vs_t()
    with open('{}/{}'.format(cur_dir, 'msd_vs_t-test0.csv'), 'w') as wf:
        for key in sorted(data.keys()):
            dr2, da2, dx, dy, da, dxdy, dxda, dyda, count = data[key][0], data[key][1], data[key][
                2], data[key][3], data[key][4], data[key][5], data[key][6], data[key][7], data[key][8]
            wf.write('{},{},{},{},{},{},{},{},{}'.format(key * dt, dr2 / count, da2 / count,
                                                         dx / count, dy / count, da / count, dxdy / count, dxda / count, dyda / count))
            wf.write('\n')
