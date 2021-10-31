import os
import sys
# import random
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

cur_dir = os.getcwd()

test1, test2 = [], []


def load(file_name):
  data = np.empty((0, 9))
  with open('{}/{}'.format(cur_dir, file_name), 'r') as rf:
    for line in rf:
      line = line.split(',')
      t, dr2, da2, dx, dy, da, dxdy, dxda, dyda = float(line[0]), float(line[1]), float(line[2]), float(
          line[3]), float(line[4]), float(line[5]), float(line[6]), float(line[7]), float(line[8])
      dr2 -= (dx * dx + dy * dy)
      da2 -= da * da
      dxda -= dx * da * 1
      dyda -= dy * da * 1
      # if dx * da < 0:
      #   print(dx, da)
      # print(dxda)
      test2.append(dxda)
      # print('------------')

      data = np.vstack(
          (data, np.array([t, dr2, da2, dx, dy, da, dxdy, dxda, dyda])))
  return data


# data = load('msd_vs_t-test0.csv')
fig = plt.figure(figsize=(15, 7), dpi=100)
gs = gridspec.GridSpec(2, 3, wspace=0.25, hspace=0.15,
                       top=0.99, bottom=0.10, left=0.05, right=0.99)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])
ax4 = fig.add_subplot(gs[3])
ax5 = fig.add_subplot(gs[4])
ax6 = fig.add_subplot(gs[5])

trials = [_ for _ in os.listdir(cur_dir) if 'msd_vs_t-test0' in _]


ax1.plot(np.linspace(0, 30, 30), 0.06 *
         np.linspace(0, 30, 30)**(1), linestyle='--', label=r'$t^{1}$', c='k', linewidth=2)

ax1.plot(np.linspace(10, 1000, 100), 0.3 *
         np.linspace(10, 1000, 100)**(0.5), linestyle='-.', label=r'$t^{0.5}$', c='r', linewidth=2)

ax2.plot(np.linspace(0, 1000, 100), 1 *
         np.linspace(0, 1000, 100)**(1), linestyle='--', label=r'$t^{1}$', c='k', linewidth=2)

ax1.text(0.60, 0.15,
         r'<drdr>($m^2$)', verticalalignment='top', horizontalalignment='left', transform=ax1.transAxes, fontsize=12, style='italic')
ax2.text(0.05, 0.90,
         r'<d$\theta$d$\theta$>($rad^2$)', verticalalignment='top', horizontalalignment='left', transform=ax2.transAxes, fontsize=12, style='italic')
ax3.text(0.05, 0.90,
         r'<dx>($m$) & <dy>($m$)', verticalalignment='top', horizontalalignment='left', transform=ax3.transAxes, fontsize=12, style='italic')
ax4.text(0.05, 0.15,
         r'<d$\theta$>($rad$)', verticalalignment='top', horizontalalignment='left', transform=ax4.transAxes, fontsize=12, style='italic')
ax5.text(0.60, 0.15,
         r'<dxdy>($m^2$)', verticalalignment='top', horizontalalignment='left', transform=ax5.transAxes, fontsize=12, style='italic')
ax6.text(0.05, 0.15,
         r'<dxd$\theta$>($m rad$)', verticalalignment='top', horizontalalignment='left', transform=ax6.transAxes, fontsize=12, style='italic')


for trial in trials:
  data = load(trial)
  time = data[:, 0]
  print(time)

  ax1.plot(time, data[:, 1], label=trial[len(trial) - 8:len(trial) - 4])
  ax1.set_xscale('log')
  ax1.set_yscale('log')

  ax2.plot(time, data[:, 2], label=r'<d$\theta$d$\theta$>($rad^2$)')
  ax2.set_xscale('log')
  ax2.set_yscale('log')

  ax3.plot(time, data[:, 3], label=r'<dx>($m$)')
  ax3.plot(time, data[:, 4], label=r'<dy>($m$)')

  ax4.plot(time, data[:, 5], label=r'<d$\theta$>($rad$)')
  ax4.set_xlabel('t(s)', fontdict={'fontsize': 20})

  ax5.plot(time, data[:, 6], label=trial[len(trial) - 8:len(trial) - 4])
  ax5.set_xlabel('t(s)', fontdict={'fontsize': 20})

  ax6.plot(time, data[:, 7], label=trial[len(trial) - 8:len(trial) - 4])
  ax6.set_xlabel('t(s)', fontdict={'fontsize': 20})

  # ax6.plot(np.linspace(0.5, 30, 30), 0.003 *
  #          np.linspace(0.5, 30, 30)**(1), linestyle='--', label=r'$t^{1}$', c='k', linewidth=2)

  # ax6.plot(np.linspace(50, 1000, 100), 0.0001 *
  #          np.linspace(50, 1000, 100)**(1.8), linestyle='-.', label=r'$t^{1.8}$', c='r', linewidth=2)

  ax1.tick_params(labelsize=15)
  ax1.legend(fontsize=10)
  ax2.tick_params(labelsize=15)
  # ax2.legend(fontsize=10)
  ax3.tick_params(labelsize=15)
  # ax3.legend(fontsize=10)
  ax4.tick_params(labelsize=15)
  # ax4.legend(fontsize=10)
  ax5.tick_params(labelsize=15)
  ax5.legend(fontsize=10)
  ax6.tick_params(labelsize=15)
  ax6.legend(fontsize=10)
# ax2.set_xlim([10**0, 10**2])
# ax2.set_ylim([0.8 * 10**-4, 10**2])
plt.show()
# fig.savefig('./fig3-2.png')
