import sys 
import matplotlib.pyplot as plt
import matplotlib
import numpy as num


font = {'family' : 'normal',
        'size'   : 9}
matplotlib.rc('font', **font)

use_scatter = True
scatter_type = 'depth_location'
use_abs = True
use_depth_difference = True
only_failed = False
xlabel = 'Mislocalization [km]'
ylabel = 'Misangle [deg]'
suptitle = ''
correct_depth = 2000
#correct_depth = 5000
print 'CORRECT DEPTH_________________ ', correct_depth
grace = 200
cmap = matplotlib.cm.get_cmap('jet')

if correct_depth==2000:
    vmin = -3
    vmax = 3
if correct_depth==5000:
    vmin = -5
    vmax = 5
if scatter_type=='depth_location':
    vmin=None
    vmax=None

print 'VMIN VMAX____________________', vmin, vmax

file_name = sys.argv[1]

f = open(file_name, 'r')
results = []

for l in f.readlines():
    results.append(map(float, l.split()))

fig = plt.figure(figsize=(4,3), dpi=100) #, frameon=False, tight_layout=True)
ax = fig.add_subplot(111)
max_data = 790 
print 'max data: ', max_data
i=1

zmin = min(num.array(results).T[3])
zmax = max(num.array(results).T[3])


cnorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
#cnorm = matplotlib.colors.Normalize(vmin=zmin-0.1*zmin, vmax=zmax+0.1*zmax)
scalarMap = matplotlib.cm.ScalarMappable(norm=cnorm, cmap=cmap)

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
        ax.plot(a,b, c, markersize=3.3)
    elif d==0.:
        c = 'ro'
        ax.plot(a,b, c, markersize=3.3)
    else:
        use_scatter = True
        break
    i+=1

print 'total got it: ', gotit
print '%s percent got it '%(float(gotit)/float(i)*100)

results = results[:max_data]

if only_failed:
    print 'only failed'
    results_gotit = []
    results_no = []
    for r in results:
        if r[3]==1 or abs(r[3]-correct_depth)<=grace:
            results_gotit.append(r)
        else:
            results_no.append(r)
    results = num.array(results_no)
    results_gotit = num.array(results_gotit)
    if use_abs:
        sc = ax.scatter(abs(results_gotit.T[0]), abs(results_gotit.T[1]),
            c='0.75', s=8, lw=0.5, alpha=0.5)
    else:
        sc = ax.scatter(results_gotit.T[0], results_gotit.T[1],
            c='0.75', s=8, lw=0.5, alpha=0.5)

else:
    results = num.array(results)

if use_scatter:
    print 'scatterplot'
    print 'scatter type: ', scatter_type

    if scatter_type=='angle_location':
        X = results.T[0]
        Y = results.T[1]
        Z = results.T[3]
        if use_abs:
            X = abs(X)
            Y = abs(Y) 
        if use_depth_difference:
            Z-=correct_depth
            cb_label = 'depth difference [km]'
        else:
            cb_label = 'z [km]'
        Z/=1000
        Z*=-1
        X/=1000

    elif scatter_type=='depth_location':
        X = results.T[0]
        Y = results.T[3]
        Z = results.T[1]
        if use_abs:
            X = abs(X)
            Z = abs(Z) 

        vmin = min(Z) 
        vmax = max(Z) 
        cb_label = 'angle [deg]'
        
    sc = ax.scatter(X, Y, c=Z, s=8, lw=0.2, vmin=vmin, vmax=vmax, cmap=cmap)
plt.colorbar(sc, label=cb_label)
    

print 'total number of tests: ', i
typestr = ''
if use_scatter:
    typestr+='_zccode'
if only_failed:
    typestr+= '_only_failed'

if scatter_type == 'depth_location':
    plt.ylim([0, 50])
    plt.xlim([0, 15])

plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.suptitle(suptitle)
plt.savefig('%s%s.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')

histfig = plt.figure(figsize=(4,3), dpi=100)
hax = histfig.add_subplot(111)
if only_failed:
    concat = num.concatenate((results_gotit, results_no))
else:
    concat = results
depths = set(concat.T[3])
print 'depth: ', depths
hax.hist(concat.T[3],len(depths))

plt.savefig('%s%s_his.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')

plt.show()
