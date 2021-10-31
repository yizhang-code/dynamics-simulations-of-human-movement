import os
import sys

cur_dir = os.getcwd()
subdirs = [_ for _ in os.listdir(cur_dir) if 'step' in _]


out = {}
for subdir in subdirs:
    T = subdir[subdir.find('-T') + 2:subdir.find('-unit')]
    with open('{}/{}/output.csv'.format(cur_dir, subdir), 'r') as rf:
        for line in rf:
            pass
        line = T + ',' + line
        out[float(T)] = line

with open('{}/cur-output.csv'.format(cur_dir), 'w') as wf:
    for T, line in sorted(out.items()):
        wf.write(line)
