#/usr/bin/python2.7
import sys 
import matplotlib.pyplot as plt
import matplotlib
import numpy as num
from scipy.interpolate import griddata, Rbf
from scipy import signal
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.ticker import FuncFormatter

font = {'family' : 'normal',
        'size'   : 9}
matplotlib.rc('font', **font)

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] == True:
        return s + r'$\%$'
    else:
        return s + '%'

def projected_2quad(points):
    plen = len(points)
    A = num.ones((plen*4, 2))
    A[:][:] = num.NAN
    # upper right
    A[:plen] = points
    # lower left
    A[plen:plen*2] = -points
    # bottom right 
    points.T[0]*=-1
    A[plen*2:plen*3] = points
    # lower right 
    points *=-1
    A[plen*3:plen*4] = points

    return A


def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = num.mgrid[-size:size+1, -sizey:sizey+1]
    g = num.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()


def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)

def plot_point_cov(points, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).

    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse(cov, pos, nstd, ax, **kwargs)

def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip


use_scatter = True
scatter_type = 'angle_location'
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


if correct_depth==2000:
    vmin = -3
    vmax = 3
if correct_depth==5000:
    vmin = -5
    vmax = 5
if scatter_type=='depth_location':
    vmin=None
    vmax=None
dz=0.2
#cmap = matplotlib.cm.get_cmap('jet')
#cmap = matplotlib.cm.get_cmap('brg')
#clrs = 'rgb'
rgba01 = ((1,0,0), (0,1,0), (0,0,1))
c_converter = matplotlib.colors.ColorConverter()
clrs = c_converter.to_rgba_array(rgba01, alpha=0.75)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(colors=clrs, 
                                                           name='my_cmap', 
                                                           gamma=1.0)
print 'VMIN VMAX____________________', vmin, vmax

file_name = sys.argv[1]

f = open(file_name, 'r')
results = []

for l in f.readlines():
    results.append(map(float, l.split()))

fig = plt.figure(figsize=(7,3), dpi=100) #, frameon=False, tight_layout=True)
ax = fig.add_subplot(121)
max_data = 1000
print 'max data: ', max_data
i=1

zmin = min(num.array(results).T[3])
zmax = max(num.array(results).T[3])


cnorm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
#cnorm = matplotlib.colors.Normalize(vmin=zmin-0.1*zmin, vmax=zmax+0.1*zmax)
scalarMap = matplotlib.cm.ScalarMappable(norm=cnorm, cmap=cmap)

cb = None
gotit = 0

#for a,b,mf,d in results[:max_data]:
#   if use_abs:
#       a = abs(a)
#       b = abs(b)
#
#   if not d in [0,1]:
#       if not use_scatter or not isinstance(d, float):
#           if abs(d-correct_depth)<=grace:
#               d=1
#           else:
#               d=0
#   
#   if d==1.:
#       gotit+=1
#       c = 'bo'
#       ax.plot(a,b, c, markersize=3.3)
#   elif d==0.:
#       c = 'ro'
#       ax.plot(a,b, c, markersize=3.3)
#   else:
#       use_scatter = True
#       break
#   i+=1

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

        #fine_points = num.where(num.abs(results.T[0]-correct_depth)<= grace)



        try:
            scaling = results.T[4]
        except IndexError:
            scaling = None

        if use_abs:
            X = num.abs(X)
            Y = num.abs(Y) 
        if use_depth_difference:
            Z-=correct_depth
            cb_label = 'vertical upshift [km]'
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
            X = num.abs(X)
            Z = num.abs(Z) 

        vmin = Z.min()
        vmax = Z.max()
        cb_label = 'angle [deg]'
        
    sc = ax.scatter(X, Y, c=Z, s=8, lw=0.2, vmin=vmin, vmax=vmax, cmap=cmap)
    projected = projected_2quad(num.array([X,Y]).T)
    plot_point_cov(projected, nstd=1, alpha=0.2, facecolor='green', 
                  edgecolor='black',
                  zorder=0)

bounds = num.arange(vmin, vmax, dz)
plt.colorbar(sc, label=cb_label, boundaries=bounds)


print 'total number of tests: ', i
typestr = ''
if use_scatter:
    typestr+='_zccode'
if only_failed:
    typestr+= '_only_failed'

if scatter_type == 'angle_location':
    plt.ylim([0, 100])
    plt.xlim([0, X.max()])
    #print 'DEACTIVATED X/Y LIMS for TESTING!'

plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.suptitle(suptitle)
plt.savefig('%s%s.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')

if scaling is not None:
    figscaling = plt.figure(figsize=(4,3), dpi=100) #, frameon=False, tight_layout=True)
    axscaling = figscaling.add_subplot(111)
    
    # Grid data:
    xg, yg = num.mgrid[X.min():X.max():100j, Y.min():Y.max():100j]

    vg = griddata((X,Y), scaling, (xg,yg), method='cubic')
    #blur_image(vg,3)
    #rbf = Rbf(X,Y,Z, epsilon=2)
    #vg = rbf(xg, yg)


    sc = axscaling.scatter(X,Y, c=scaling, s=8, lw=0.2, zorder=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar(sc, label='scaling factor')
    #plt.colorbar(sc, label='misfit M')
    plt.savefig('%s%s_scaling.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')
    plt.contourf(xg,yg,vg, zorder=0)

    
plt.ylim([0, 100])
plt.xlim([0, X.max()])

#histfig = plt.figure(figsize=(4,3), dpi=100)
hax = fig.add_subplot(122)
if only_failed:
    concat = num.concatenate((results_gotit, results_no))
else:
    concat = results
depths = set(concat.T[3])
print 'depth: ', depths
hax.hist(concat.T[3],
         len(depths), 
         orientation='horizontal',
         #range=[bounds.min(), bounds.max()],
         facecolor='gray', 
         normed=1,
         align='mid')

formatter = FuncFormatter(to_percent)

plt.gca().xaxis.set_major_formatter(formatter)
plt.gcf().autofmt_xdate()

plt.xlabel('vertical mislocation [km]')
plt.ylabel('number')

plt.savefig('%s%s_his.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')

plt.show()
