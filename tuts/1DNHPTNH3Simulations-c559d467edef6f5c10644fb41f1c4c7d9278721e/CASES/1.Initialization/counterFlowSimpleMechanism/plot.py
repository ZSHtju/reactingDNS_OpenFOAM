# !/usr/bin/python
# coding: utf8
# @Time    : 2022-07-03 21:36:18
# @Author  : Shijie Xu
# @Email   : shijie.xu@energy.lth.se
# @Software: OpenFOAM

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#import pyFoam as pf

plt.style.use(['science','no-latex','ieee'])

import matplotlib as mpl
mpl.rcParams['axes.facecolor'] = 'none'
mpl.rcParams['figure.facecolor'] = 'none'
#mpl.rcParams['font.family']=['serif']
#mpl.rcParams['font.serif']=['Times New Roman']
#mpl.rcParams['lines.linewidth']=1.0
mpl.rcParams['lines.markersize']=3.0
mpl.rcParams['markers.fillstyle']='none'
#mpl.rcParams['figure.figsize'] = [3.3, 2.5]
#mpl.rcParams['font.size'] = 8.0
#mpl.rcParams['legend.fontsize'] = 'medium'
#mpl.rcParams['figure.titlesize'] = 'large'
#mpl.rcParams['figure.dpi'] = 80
#mpl.rcParams['savefig.dpi'] = 600
#ax.ticklabel_format(style='sci', scilimits=(-1,2), axis='y')
fig, ax1 = plt.subplots()

EXP=pd.read_csv('postProcessing/sampleSets/0.75/data_NH3_C7H16_O2_N2_p_T.csv')

color = 'tab:red'
ax1.tick_params(axis='y', labelcolor=color)
ax1.plot(EXP['x']*1000, EXP['T'], color='tab:red', label='T')
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.tick_params(axis='y')


ax2.plot(EXP['x']*1000, EXP['NH3'], label='NH3')
ax2.plot(EXP['x']*1000, EXP['C7H16'], label='C7H16')
ax2.plot(EXP['x']*1000, EXP['O2'], label='O2')
ax2.plot(EXP['x']*1000, EXP['N2'], label='N2')
plt.legend(loc='best')

ax1.set_xlabel("x [mm]")
ax1.set_ylabel("Temperature [K]")
ax2.set_ylabel("Species mass fractions [-]")
plt.tick_params(direction='in')
plt.tight_layout()

plt.savefig('mixingLayerProfiles.pdf', format='pdf')

