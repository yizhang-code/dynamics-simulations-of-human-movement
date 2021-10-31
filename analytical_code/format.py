import numpy as np
import math
import os
import sys
import csv
import concurrent.futures

cur_dir = os.getcwd()


def reduce_filesize(kid):
    with open(cur_dir + '/info_' + str(kid) + '.csv', 'r') as rf:
        with open(cur_dir + '/xy_' + str(kid) + '.csv', 'w') as wf:
            csv_writer = csv.writer(wf, delimiter=',')
            for i, line in enumerate(rf):
                # print(i, line)
                if i % chunk_unit_size == 0:
                    line = line.strip().split(',')
                    # print(line)
                    # print(float(line[0]), float(line[1]), float(line[2]))
                    cx = float(line[0])
                    cy = float(line[1])
                    agl = float(line[2])

                    # print(cx, cy, agl)
                    lx, ly = np.array([cx, cy]) + 0.2 * \
                        np.array([-math.sin(agl), math.cos(agl)])
                    rx, ry = np.array([cx, cy]) + 0.2 * \
                        np.array([math.sin(agl), -math.cos(agl)])
                    csv_writer.writerow(
                        (i // chunk_unit_size, lx, ly, rx, ry))
                # output(kid, i // chunk_unit_size, lx, ly, rx, ry)

    return True


kids = len([_ for _ in os.listdir(cur_dir) if 'info_' in _])
chunk_unit_size = int(sys.argv[1])
if __name__ == '__main__':
    with concurrent.futures.ProcessPoolExecutor() as executor:
        recurrence = np.arange(0, kids, 1)
        # print(recurrence)
        results = executor.map(reduce_filesize, recurrence)
