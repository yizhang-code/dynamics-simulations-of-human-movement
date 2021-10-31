import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import os
import itertools

cur_dir = os.getcwd()

files = [_ for _ in os.listdir(cur_dir) if '.dat' in _]
# print(files)

data = {}
for file in files:
    x, y = [], []
    # print(cur_dir + '/file')
    with open(cur_dir + '/' + file, 'r') as rf:
        for _, line in enumerate(rf):
            line = line.strip().split('\t')
            print(line)
            r, gr = float(line[0]), float(line[1])
            x.append(r)
            y.append(gr)
        data[file[:len(file) - 6]] = (x, y)

# print(data)
f, ax = plt.subplots(figsize=(9.5, 7.5))
for ele in sorted(data.keys()):
    x, y = data[ele][0], data[ele][1]
    # print(x, y)
    ax.plot(x, y, label=str(ele), linestyle='-')
plt.xlabel('r/meters', fontsize=30)
plt.ylabel('g(r)', fontsize=30)
plt.axhline(y=1, linestyle='-', color='k')
plt.xticks(fontsize=25, rotation=0)
plt.yticks(fontsize=20, rotation=0)
plt.title('MC simulation', fontsize=30)
plt.legend(fontsize=15, ncol=1)
plt.xlim([-0.2, 6.0])
plt.ylim([-0.2, 4.0])
plt.show()
