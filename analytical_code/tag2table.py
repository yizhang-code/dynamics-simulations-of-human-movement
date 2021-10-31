import numpy as np
import math
import os
import sys
import csv
# import concurrent.futures

cur_dir = os.getcwd()


def output(data_matrix):
    with open(cur_dir + '/all_xy.csv', 'a') as af:
        csv_writer = csv.writer(af, delimiter=',')
        for row in data_matrix:
            csv_writer.writerow(row)


def load_file(start, end):
    data_matrix = np.zeros((end - start, kids * 4 + 1))
    for kid in range(kids):
        with open(cur_dir + '/xy_' + str(kid) + '.csv', 'r') as rf:
            for step, line in enumerate(rf):
                if step >= start and step < end:
                    line = line.split(',')
                    if kid == 0:
                        data_matrix[step - start, 0] = step
                    for i in range(4):
                        data_matrix[step - start, kid *
                                    4 + i + 1] = line[i + 1]
                pass
    return data_matrix


def initial_output():
    with open(cur_dir + '/all_xy.csv', 'w') as wf:
        pass

kids = len([_ for _ in os.listdir(cur_dir) if 'info_' in _])
start = 0
chunk_unit_size = 500
end = start + chunk_unit_size
if __name__ == '__main__':
    initial_output()
    while True:
        chunk_matrix = load_file(start, end)
        # print(mtx)
        if sum(chunk_matrix[0]) == 0:
            break
        output(chunk_matrix)
        start = end
        end += chunk_unit_size
