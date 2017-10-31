#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 12:20:07 2017

@author: zonni
"""

import matplotlib.pyplot as plt
from matplotlib import rcParams

times = [14.0, 22.5, 11.9, 9.1, 44.8, 36.8, 14.0, 16.7, 85.0,  52.2, 19.3, 14.0]
i = [241, 3771, 2129, 1509, 7447, 4542, 2484, 2914, 14264, 7876, 3364, 2484]

times.sort()
i.sort()

params = {'backend': 'pdf',
          'font.family': 'serif',
          'font.size': 16,
          'axes.labelsize': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': True,
          'lines.linewidth': 2,
          'lines.markersize': 15,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1
          }
          
rcParams.update(params)
plt.close('all')


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(i,times, '.', c='#1f77b4')

ax1.grid()
ax1.set_ylabel("CPU Time (s)")
ax1.set_xlabel("Inner Iterations")
plt.savefig("InnerDependence.pdf", format='pdf')