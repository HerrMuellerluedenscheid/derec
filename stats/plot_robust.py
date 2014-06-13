import sys 
import matplotlib.pyplot as plt
import matplotlib
import numpy as num

#cmap = matplotlib.cm.jet(20)
cmap = matplotlib.cm.get_cmap('jet')


def SaveFigureAsImage(fileName,fig=None,**kwargs):
    ''' Save a Matplotlib figure as an image without borders or frames.
       Args:
            fileName (str): String that ends in .png etc.

            fig (Matplotlib figure instance): figure you want to save as the image
        Keyword Args:
            orig_size (tuple): width, height of the original image used to maintain 
            aspect ratio.
    '''
    fig_size = fig.get_size_inches()
    w,h = fig_size[0], fig_size[1]
    fig.patch.set_alpha(0)
    if kwargs.has_key('orig_size'): # Aspect ratio scaling if required
        w,h = kwargs['orig_size']
        w2,h2 = fig_size[0],fig_size[1]
        fig.set_size_inches([(w2/w)*w,(w2/w)*h])
        fig.set_dpi((w2/w)*fig.get_dpi())
    a=fig.gca()
    a.set_frame_on(False)
    a.set_xticks([]); a.set_yticks([])
    plt.axis('off')
    plt.xlim(0,h); plt.ylim(w,0)
    fig.savefig(fileName, transparent=True, bbox_inches='tight', \
                        pad_inches=0)



font = {'family' : 'normal',
        'size'   : 9}

use_scatter = True 
use_abs = False 
only_failed = True 

matplotlib.rc('font', **font)

file_name = sys.argv[1]

f = open(file_name, 'r')
results = []

for l in f.readlines():
    results.append(map(float, l.split()))

# for plotting 1 of 2 figures side by side on a4:

fig = plt.figure(figsize=(4,3), dpi=100) #, frameon=False, tight_layout=True)
#fig = plt.figure()
ax = fig.add_subplot(111)
max_data = 2000
i=1

zmin = min(num.array(results).T[3])
zmax = max(num.array(results).T[3])
print zmin, zmax


cnorm = matplotlib.colors.Normalize(vmin=zmin-0.1*zmin, vmax=zmax+0.1*zmax)
scalarMap = matplotlib.cm.ScalarMappable(norm=cnorm, cmap=cmap)
correct_depth = 5000
#correct_depth = None
grace = 200

cb = None
gotit = 0

for a,b,mf,d in results[:max_data]:
    if use_abs:
        a = abs(a)
        b = abs(b)
        
    if not d in [0,1]:
        if not use_scatter or not isinstance(d, float):
            if abs(d-correct_depth)<=grace:
                d=1
            else:
                d=0
    
    if d==1.:
        gotit+=1
        c = 'bo'
        ax.plot(a/1000.,b, c, markersize=3.3)
    elif d==0.:
        c = 'ro'
        ax.plot(a/1000.,b, c, markersize=3.3)
    else:
        use_scatter = True
        break
    i+=1

print 'total got it: ', gotit
print '%s percent got it '%(float(gotit)/float(i)*100)

results = results[:max_data]

if only_failed:
    results_gotit = []
    results_no = []
    for r in results:
        if r[3]==1 or abs(r[3]-correct_depth)<=grace:
            results_gotit.append(r)
        else:
            results_no.append(r)
    results = num.array(results_no)
    results_gotit = num.array(results_gotit)
    sc = ax.scatter(results_gotit.T[0], results_gotit.T[1],
            c='0.75', s=8, lw=0.5, alpha=0.5)
else:
    results = num.array(results)

if use_abs:
    map(abs, results.T[0])
    map(abs, results.T[1])

    
if use_scatter:
    print 'scatterplot'

    sc = ax.scatter(results.T[0], results.T[1],
            c=results.T[3], s=8, lw=0.2,
            vmin=zmin, vmax=zmax, cmap=cmap)
    plt.colorbar(sc, label='z [m]')
    

#color = scalarMap.to_rgba(d, alpha=0.5)
#c='o'
#ax.plot(abs(a)/1000.,abs(b), c, color=color, markersize=3.3)
#if not cb:
#    i+=1

print 'total number of tests: ', i
#plt.colorbar(im)
typestr = ''
if use_scatter:
    typestr+='_zccode'
if only_failed:
    typestr+= '_only_failed'

#plt.xlabel('Mislocation [km]' )
#plt.ylabel('Mislocation [km]' )
#plt.ylabel('Angle [deg]')
#plt.xlim([0, 12])
#plt.ylim([0, 60])
#plt.suptitle(file_name)
plt.savefig('%s%s.pdf'%(file_name.split('.')[0], typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')

histfig = plt.figure(figsize=(4,3), dpi=100)
hax = histfig.add_subplot(111)
if only_failed:
    concat = num.concatenate((results_gotit, results_no))
else:
    concat = results
depths = set(concat.T[3])
print depths
hax.hist(concat.T[3], len(depths))
plt.savefig('%s%s_his.pdf'%(file_name.split('.')[0], typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')

plt.show()
