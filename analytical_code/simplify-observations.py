# import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import os

cur_dir = os.getcwd()
tar_dir = sys.argv[1]
dst_dir = '{}/{}'.format(cur_dir, tar_dir)
window_size = 5


class Smooth1Dlist(object):
    """docstring for SmoothData"""

    def __init__(self, data):
        self.data = data
        self.window_size = window_size

    def extend_list(self, ori_list):
        extend_head, extend_tail = [], []
        i = self.window_size
        while i > 0:
            extend_head.append(ori_list[i])
            extend_tail.append(ori_list[len(ori_list) - window_size + i - 2])
            i -= 1
        return extend_head + ori_list + extend_tail

    @staticmethod
    def ave_items(sub_list):
        tot_val, tot_count = 0, 0
        for item in sub_list:
            val, count = item
            tot_val += val
            tot_count += count
        # print(tot_val, tot_count)
        return tot_val / tot_count

    def smooth(self):
        cos_list = []
        sorted_r = sorted(self.data.keys())
        for r in sorted_r:
            cos_list.append(self.data[r])
        ori_data_length = len(cos_list)
        cos_list = self.extend_list(cos_list)
        smooth_cos_list = []
        for i in range(ori_data_length):
            smooth_cos_list.append(self.ave_items(
                cos_list[i:i + 2 * window_size + 1]))
        return sorted_r, smooth_cos_list


class LoadFiles(object):
    """docstring for LoadFiles"""

    def __init__(self, observation):
        self.observation = observation
        self.observation_dict = {}

    def add_item(self, line):
        r, agl1, agl2 = float(line[0]), float(line[1]), float(line[2])
        r = round(int(r / 0.08) * 0.08, 2)
        cos_val = - math.cos((agl1 - agl2) * np.pi / 180)
        if r not in self.observation_dict:
            self.observation_dict[r] = [cos_val, 1]
        else:
            self.observation_dict[r][0] += cos_val
            self.observation_dict[r][1] += 1

    def read_files(self):
        # for file in self.files:
        with open('{}/{}'.format(dst_dir, observation), 'r') as rf:
            for _, line in enumerate(rf):
                line = line.split('\t')
                self.add_item(line)

        smooth_r, smooth_cos = Smooth1Dlist(self.observation_dict).smooth()
        with open('{}/smooth_{}.dat'.format(cur_dir, self.observation[:len(self.observation) - 4]), 'w') as wf:
            print('start to write {} info to file'.format(
                self.observation[:len(self.observation) - 4]) + '\n')
            for r, cos_val in zip(smooth_r, smooth_cos):
                wf.write(str(r) + '\t' + str(cos_val))
                wf.write('\n')
        print('{} info has written to file!'.format(
            self.observation[:len(self.observation) - 4]) + '\n')
        return True


if __name__ == '__main__':
    observations = [_ for _ in os.listdir(dst_dir) if 'dat' in _]
    for observation in observations:
        data_obj = LoadFiles(observation)
        data_obj.read_files()
