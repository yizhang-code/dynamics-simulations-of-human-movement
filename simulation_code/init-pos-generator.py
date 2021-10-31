import random
from math import sqrt
import numpy as np
import os
import sys
from shutil import copy

cur_dir = os.getcwd()

init_length = 8.5
init_width = 4.8

num = int(sys.argv[2])

density = float(sys.argv[1])

radius = 0.25
length = np.sqrt(np.pi * radius * radius * num / density)
width = length


unit_num = int(np.sqrt(num)) + 1


def random_generator(length, width, num):
    loc_list = []
    for i in range(0, unit_num):
        for j in range(0, unit_num):
            xya = (length * i / unit_num, width * j / unit_num, 0)
            loc_list.append(xya)
            if len(loc_list) >= num:
                return loc_list
    if len(loc_list) < num:
        print('Size is too small!!!')
        return []


loc_list = random_generator(length, width, num)
print(length, width, num)
print('Inital density ------------>{}\nNew density --------------->{}'.format(10 /
                                                                              8.5 / 4.8, np.pi * num * radius * radius / length / width))


with open('{}/initial_position-many.dat'.format(os.getcwd()), 'w') as wf:
    for item in sorted(loc_list):
        wf.write('{}\t{}\t{}'.format(item[0], item[1], item[2]))
        wf.write('\n')

try:
    print(int(density))
    os.mkdir('{}/density-{}'.format(cur_dir, int(density)))
except:
    pass
try:
    copy('{}/initial_position-many.dat'.format(cur_dir),
         '{}/density-{}/initial_position-many.dat'.format(cur_dir, int(density)))
except:
    print('Err ---- can not copy position file!!')
