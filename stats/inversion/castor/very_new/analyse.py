import os
import sys
import shutil
import glob
from pyrocko import mopad,moment_tensor
from collections import defaultdict
import pylab
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as num

def get_angles(mts):
    out = []
    omt = mts[0].get_pyrocko_moment_tensor()
    for mt in mts:
        try:
            a = mt.get_pyrocko_moment_tensor().angle(omt)
        except:
            a = 0.
        out.append(a)
    return out

max_runs = 15
fns = sys.argv[1:]
numsteps = 0
fsize= 8
results = defaultdict()
print fns
for fn in fns:
    f = open(fn, 'r')
    tmp = []
    for l in f.readlines():
        tmp.append(l.split())
    print tmp

    results[fn] = tmp

dirname = 'tmp'
if not os.path.exists(dirname):
    os.mkdir(dirname)

elapsed = defaultdict()
mfs = defaultdict()
fnouts = []
angles = defaultdict()
irun=0
for fn,result in results.items():
    if irun==max_runs:
        break
    print 'ASDFASDF', irun
    localdirname = fn.split('.')[0]
    if not os.path.exists(localdirname):
        os.mkdir(localdirname)
    mf = []
    mts = []
    numsteps = max(numsteps, len(result)-1)
    for i,step in enumerate(result[1::]):


        print step
        mf.append(float(step[0]))
        #z =float(step[1])
        s =float(step[1])
        d =float(step[2])
        r =float(step[3])
        
        if s<0. or s>360 or d< -90 or d>90 or r >180 or r<-180:
            results[fn]=None
            localdirname = 'delme'
            shutil.rmtree(fn.split('.')[0])
            break


        print d
        mt = mopad.MomentTensor(M=(s,d,r))
        mts.append(mt)

        bb = mopad.BeachBall(mt)
        f = bb.ploBB(kwargs={'_return_fig':True})
        print f
        fnout = fn.split('.')[0]+'_step%s.png'%str(i)
        f.set_dpi(400)
        f.savefig(localdirname+'/'+fnout, pad_inches=0.0, frameon=False,
                bbox_inches='tight')
    
    if localdirname!='delme':
        angles[localdirname] = get_angles(mts)
        mfs[localdirname] = mf
        elapsed[localdirname] = float(result[0][2])
        fnouts.append(localdirname)
    irun+=1

f, axs = plt.subplots(numsteps, len(fnouts)+1 , dpi=400)
figwith = 6.
figheight = numsteps/(len(fnouts)+1.)*figwith
f.set_figwidth(figwith)
f.set_figheight(figheight)
boxprops = dict(boxstyle='round', facecolor='white', alpha=0.5, linewidth=0.)

for fni, fn in enumerate(fnouts[:max_runs:]):
    ela = elapsed[fn]
    imgfns = glob.glob(fn+'/*')
    imgs = defaultdict()
    for imgfn in imgfns:
        step = imgfn.split('__')[1]
        step = int(step[4])
        img = mpimg.imread(imgfn)
        imgs[step] = img
    
    for i, img in imgs.items():
        try:
            ax = axs[i-1, fni+1]
        except IndexError:
            break 
        ax.imshow(img, interpolation='bilinear')

        #ax.text(0,1, '%1.2f'%(angles[fn][i]), horizontalalignment='left',
        #        verticalalignment='top', transform=ax.transAxes, fontsize=fsize)
        if not mfs[fn][i]==0:
            ax.text(1.,1., '%1.2f'%(mfs[fn][i]), horizontalalignment='right',
                    verticalalignment='top', transform=ax.transAxes,
                    fontsize=fsize-1,
                    bbox=boxprops)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.axis('off')

toprow = axs[:,0]
for i in range(numsteps-1):
    toprow[i].text(-0.1,0.5,'Step %s'%(i+1), fontsize=fsize)
toprow[numsteps-1].text(-0.1,0.5,'Input', fontsize=fsize)
toprow[numsteps-1].text(0.,0, 'mean elapsed time: %1.2f s'%(num.mean(ela)), horizontalalignment='left',
    verticalalignment='top', transform=toprow[numsteps-1].transAxes,
    fontsize=fsize)
for ax in toprow:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.patch.set_visible(False)
    ax.axis('off')

f.subplots_adjust(hspace=0.0, wspace=0.1)
f.savefig('bb_compiled.pdf', transparent=False, pad_inches=0.10, dpi=400,
        bbox_inches='tight')
                        
plt.show()
 

    



