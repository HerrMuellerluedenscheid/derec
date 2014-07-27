#/usr/bin/python2.7
import sys 
import os
import matplotlib.pyplot as plt
import matplotlib
import numpy as num
from scipy import signal
import numpy as np

import matplotlib
matplotlib.use('PDF')
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.ticker import FuncFormatter
from matplotlib.font_manager import FontProperties
from matplotlib import rc
from collections import defaultdict
from plot_utils import *

font = {'family' : 'normal',
        'size'   : 8}
matplotlib.rc('font', **font)
plt.rcParams['ps.useafm'] = True

file_names = sys.argv[1:]
cwd = os.getcwd()
if iscastor_dir(cwd):
    correct_depth = 2.
    vmin = -3
    vmax = 3

if isdoctar_dir(cwd):
    correct_depth = 5.
    vmin = -2
    vmax = 5
        
print 'CORRECT DEPTH_________________ ', correct_depth
if isplain_dir(cwd):
    graces = [0.2]
else:
    graces = [0.2, 0.4, 0.6]
print 'graces ', graces

dz=0.2
rgba01 = ((1,0,0), (0,1,0), (0,0,1))
c_converter = matplotlib.colors.ColorConverter()
clrs = c_converter.to_rgba_array(rgba01, alpha=0.75)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(colors=clrs, 
                                                           name='my_cmap', 
                                                           gamma=1.0)
plt.register_cmap(cmap=cmap)
print 'VMIN VMAX____________________', vmin, vmax


noise_levels = [] 

#num_gotit = [] 
num_gotit = defaultdict(list)
z_stds = [] 
average_z = []

for fn in file_names:
    f = open(fn, 'r')

    results = []
    for l in f.readlines():
        results.append(map(float, l.split()))
    results = num.array(results)
    z, snr, sig = results.T[3:6]
    z /= 1000.
    map(num.array, [z,snr,sig])
    for grace in graces:
        gotit = float(len(num.where(num.abs(z-correct_depth)<=grace)[0]))
        num_gotit[grace].append(gotit*100/float(len(z)))
    z_stds.append(num.std(z))
    average_z.append(num.mean(z))
    noise_levels.append(num.mean(snr))
    sig_ag = num.mean(sig)


fig, (ax, ax2) = plt.subplots(2,1, figsize=(3,3), dpi=100, sharex=True) #, frameon=False, tight_layout=True)

#ax = fig.add_subplot(211)
#num_gotit = num.array(num_gotit)*100
#ax.plot(noise_levels, num_gotit, marker='o', linestyle='')
z_scalefactor = 0.1
z_stds = num.array(z_stds)/z_scalefactor
kwargs = {}
if len(num_gotit)==1:
    kwargs = {'c':'r'}

for grace, gotit in num_gotit.items():
    eb = ax.errorbar(noise_levels, 
                gotit, 
                marker='o', 
                linestyle='', 
                lw=1., 
                ms=4.,
                label=grace,
                **kwargs)
    
if len(num_gotit)!=1:
    ax.legend(loc='upper left',
    fontsize=8,
    numpoints=1,
    title='vert. freedom [km]', 
    bbox_to_anchor=[0,1.2], 
    ncol=len(num_gotit))

ax.set_xlim((-50,40))
ax.set_ylim((-20,110))

xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
linepos = xmax
sigma = 2.
xpos = 0.9
ypos = 0.1
xpos *=xmin+(xmax-xmin)
ypos *=ymax+(ymax-ymin)

#ax.errorbar(xpos, 
#          ypos, 
#          marker='', 
#          linestyle='', 
#          yerr=sigma/z_scalefactor,
#          c='b', 
#          lw=1.)


#ax.text(xpos, ypos, r'$\sigma_z = %skm $'%(sigma), rotation=90, ha='right', va='center')
#ax.set_xlabel('signal-noise-ratio [db]')
ax.set_ylabel('Correct depth [\%]')
yticks = ax.yaxis.get_major_ticks()
yticks[0].set_visible(False)

#fig2 = plt.figure(figsize=(4,2), dpi=100) #, frameon=False, tight_layout=True)
#ax2 = fig.add_subplot(212)
print average_z
average_z = correct_depth-num.array(average_z)
ax2.errorbar(noise_levels, average_z, marker='o', linestyle='',
        yerr=z_stds*z_scalefactor, 
        lw=0.5, zorder=1,
        ms=4.)
ax2.set_ylim([vmin, vmax])
#ax2.scatter(noise_levels, average_z, c=average_z, vmin=vmin, vmax=vmax,
#        cmap=cmap, zorder=2)
#fancy_box = dict(boxstyle='round', facecolor='wheat', alpha=1)
#ax2.axhline(y=0, linestyle='--', c='g', zorder=3, label='correct depth')
#ax2.text(x=xpos, y=correct_depth, s='correct depth',
#        verticalalignment='top',
#        horizontalalignment='right')
        #bbox=fancy_box)
ax2.set_xlabel('signal-noise-ratio [db]')
ax2.set_ylabel('Vertical upshift [km] ')
#fprop = FontProperties()
#fprop.set_size(9)
#ax2.legend(loc=0, prop=fprop)

ax.xaxis.grid(linewidth=0.5, color='grey', alpha=0.5)
ax2.xaxis.grid(linewidth=0.5, color='grey', alpha=0.5)

plt.subplots_adjust(hspace=0.1)
plt.savefig('snr_ratio.pdf', transparent=True, pad_inches=0.01, bbox_inches='tight')

plt.show()
