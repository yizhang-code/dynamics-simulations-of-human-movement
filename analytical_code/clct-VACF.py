import os
import sys
# import random
from math import sqrt
import numpy as np
# import matplotlib.pyplot as plt

cur_dir = os.getcwd()
dt = float(sys.argv[1])
# dt = 1
max_t = int(sys.argv[1])
time_unit = 1


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


class ClctVACF(object):

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

    @staticmethod
    def clct_vel(loc_list, ori_list):
        vels = np.empty((0, 3))
        for t in range(len(loc_list) - time_unit):
            x_t2, y_t2, ori_t2 = loc_list[t +
                                          time_unit][0], loc_list[t + time_unit][1], ori_list[t + time_unit]
            x_t1, y_t1, ori_t1 = loc_list[t][0], loc_list[t][1], ori_list[t]

            vt = np.array([(x_t2 - x_t1) / time_unit / dt, (y_t2 - y_t1) /
                           time_unit / dt, (ori_t2 - ori_t1) / time_unit / dt])
            vels = np.vstack((vels, vt))
        return vels

    def clctvacf(self, t):
        sum_vtv02, sum_vtv02_ori, N = 0, 0, 0
        for flag, loc_list in self.chunk_loc.items():
            ori_list = self.map2range(self.chunk_ori[flag])
            if len(ori_list) >= 2:
                vels = self.clct_vel(loc_list, ori_list)
            else:
                vels = []
            if len(vels) > t:
                for i in range(len(vels) - t):
                    dv = vels[i + t] - vels[i]
                    sum_vtv02 += vels[i + t][0] * \
                        vels[i][0] + vels[i + t][1] * vels[i][1]
                    sum_vtv02_ori = vels[i + t][2] * vels[i][2]
                    N += 1
        if N == 0:
            return (0, 0, 0)
        else:
            return (sum_vtv02 / N, sum_vtv02_ori / N, 1)

    def find_vacf(self):
        self.chunk_raw()
        for t in range(self.max_t):
            print(self.particle_num, t)
            vtv0, vtv0_ori, count = self.clctvacf(t)
            if t not in self.data:
                self.data[t] = [0, 0, 0]
            self.data[t][0] += vtv0
            self.data[t][1] += vtv0_ori
            self.data[t][2] += count
        return self.data


if __name__ == '__main__':
    data = {}
    particle_0 = LoadAllxy(0)
    particle_0_loc = particle_0.load_all_xy()
    particle_num = particle_0.particle_num
    data = ClctVACF(particle_0_loc, data, 0).find_vacf()

    for num in range(1, particle_num):
        print(num)
        data = ClctVACF(LoadAllxy(num).load_all_xy(),
                        data, num).find_vacf()
    with open('{}/{}'.format(cur_dir, 'vacf-test0.csv'), 'w') as wf:
        for key in sorted(data.keys()):
            vtv0, vtv0_ori, count = data[key][0], data[key][1], data[key][2]
            wf.write('{},{},{}'.format(
                key * dt, vtv0 / count, vtv0_ori / count))
            wf.write('\n')
