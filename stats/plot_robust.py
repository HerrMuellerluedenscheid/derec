#/usr/bin/python2.7
import derec.derec_utils as du
import os
import sys 
import matplotlib.pyplot as plt
import matplotlib
import numpy as num
try:
    from scipy.interpolate import griddata, Rbf
    nogrid = False
except:
    print 'no gridding available'
    nogrid = True
from scipy import signal
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.ticker import FuncFormatter
from matplotlib import gridspec
from plot_utils import to_percent, iscastor_dir, isdoctar_dir

font = {'family' : 'normal',
        'size'   : 9}
matplotlib.rc('font', **font)


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


def gridded_counter(ax, X,Y,Z,xstep=None, ystep=None, numx=None, numy=None, zgrace=0.2):
    xmin, xmax = ax.get_xlim()
    xmin=xmin if xmin>=0 else 0
    ymin, ymax = ax.get_ylim()
    ymin=ymin if ymin>=0 else 0
    x_low = 0
    X_centers = []
    Y_centers = []
    Z_centers = []
    if numx and numy:
        x_vals = num.linspace(xmin, xmax, numx)
        y_vals = num.linspace(ymin, ymax, numy)
        xstep = x_vals[1]
        ystep = y_vals[1]

    while x_low<xmax-xstep:
        x_up = x_low+xstep
        if x_low==xmin:
            x_center = xmin
        else:
            x_center = x_low+(x_up-x_low)/2.
        y_low = 0
        while y_low<ymax-ystep:
            y_up = y_low+ystep
            if y_low==ymin:
                y_center=ymin
            else:
                y_center = y_low+(y_up-y_low)/2.
            loc_x = num.where(num.logical_and(X<=x_up, X>=x_low))
            loc_y = num.where(num.logical_and(Y<=y_up, Y>=y_low))
            loc = num.intersect1d(loc_x[0], loc_y[0])
            if len(loc)<=2: 
                gotit_perc=0.#num.NAN
            else:
                points = Z[loc]
                gotit_perc = 1.*len(num.where(num.abs(points)<=zgrace)[0])/len(points)*100.
            
            X_centers.append(x_center)
            Y_centers.append(y_center)
            Z_centers.append(gotit_perc)
            
            y_low+=ystep
        x_low+=xstep

    XC = num.array(X_centers)
    YC = num.array(Y_centers)
    ZC = num.array(Z_centers)
    return XC, YC, ZC

cwd = os.getcwd()
use_scatter = True
scatter_type = 'angle_location'
use_abs = True
use_depth_difference = True
only_failed = False
xlabel = 'Mislocalization [km]'
ylabel = 'Misangle [deg]'
suptitle = ''
if iscastor_dir(cwd):
    correct_depth = 2000
    counterx = 0.8
    countery = 10
    x_max = 6
    y_max = 90
if isdoctar_dir(cwd):
    correct_depth = 5000
    counterx = 0.5
    countery = 10
    x_max = 4
    y_max = 90


grace = 200
graces = [0.2, 0.4]
isolevel = 66.6
if correct_depth==2000:
    vmin = -3.0
    vmax = 1.4
    #vmin = -1
    #vmax = 1
if correct_depth==5000:
    vmin = -3
    vmax = 4
if scatter_type=='depth_location':
    vmin=None
    vmax=None

gamma = num.abs(vmin/vmax)
dz=0.2
cmap_jet = matplotlib.cm.get_cmap('jet')
#cmap = matplotlib.cm.get_cmap('brg')
#clrs = 'rgb'
rgba01 = ((1,0,0), (0,1,0), (0,0,1))
c_converter = matplotlib.colors.ColorConverter()
clrs = c_converter.to_rgba_array(rgba01, alpha=0.75)
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list(colors=clrs, 
#                                                           name='my_cmap', 
#                                                           gamma=1.)

print 'VMIN VMAX____________________', vmin, vmax

file_name = sys.argv[1]
target_depth = None
if len(sys.argv)>2:
    target_depth = sys.argv[2]
    target_zmin = sys.argv[3]
    target_zmax = sys.argv[4]

if target_depth:
    correct_depth = float(target_depth)

print 'CORRECT DEPTH_________________ ', correct_depth

f = open(file_name, 'r')
results = []

for l in f.readlines():
    results.append(map(float, l.split()))

fig = plt.figure(figsize=(6,3), dpi=100, tight_layout=False) #, frameon=False, tight_layout=True)
gs = gridspec.GridSpec(1,2, width_ratios=[2,1])
ax = fig.add_subplot(gs[0])
max_data = 2000
print 'max data: ', max_data
i=1
if target_zmin:
    vmin = float(target_zmin)
if target_zmax:
    vmax = float(target_zmax )
cmap = du.get_cmap(N=len(num.array(results).T[0])+1, gamma=gamma)

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
            c='0.75', s=8, lw=0.1, alpha=0.5)
    else:
        sc = ax.scatter(results_gotit.T[0], results_gotit.T[1],
            c='0.75', s=8, lw=0.1, alpha=0.5)

else:
    results = num.array(results)

if use_scatter:
    print 'scatterplot'
    print 'scatter type: ', scatter_type

    if scatter_type=='angle_location':
        X = results.T[0]
        Y = results.T[1]
        Z = results.T[3]

        try:
            misfit = results.T[2]
        except IndexError:
            misfit = None
        try:
            scaling = results.T[4]
        except IndexError:
            scaling = None

        if use_abs:
            X = num.abs(X)
            Y = num.abs(Y) 
        if use_depth_difference:
            Z-=correct_depth
            cb_label = 'upshift [km]'
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
        
    sc = ax.scatter(X, Y, c=Z, s=8, lw=0.1, vmin=vmin, vmax=vmax, cmap=cmap, zorder=2)
    plt.ylim([0, y_max])
    plt.xlim([0, x_max])
    Xc, Yc, Zc = gridded_counter(ax, 
                                 X, 
                                 Y, 
                                 Z, 
                                 xstep=0.5, 
                                 ystep=4,  
                                 zgrace=grace/1000.)

    for zgrace in graces:
        Xc, Yc, Zc = gridded_counter(ax, X, Y, Z, xstep=counterx, ystep=countery,
                                     zgrace=zgrace)
        xg, yg = num.mgrid[Xc.min():Xc.max():50j, Yc.min():Yc.max():50j]
        if not nogrid:
            vg = griddata((Xc,Yc), Zc, (xg,yg), method='cubic')

            ax.contourf(xg,
                        yg,
                        vg,
                        linewidth=2,
                        zorder=0,
                        alpha=0.5, 
                        levels=[isolevel,110],
                        colors=('grey' )) 
            ax.contour(xg,yg,vg, levels=[isolevel], linewidths=(2), colors=('grey'), 
                      zorder=1)
    # ax.scatter(Xc, Yc, s=Zc, c='b', marker='o')

bounds = num.arange(vmin, vmax+dz, dz)
cticks = num.arange(vmin, vmax+dz, 1)
cb = plt.colorbar(sc, 
                  label=cb_label, 
                  boundaries=bounds, 
                  pad=0.1)
cb.set_ticks(cticks)
cb.set_label(cb_label, labelpad=-40)
print 'total number of tests: ', i
for t in cb.ax.get_yticklabels():
    t.set_horizontalalignment('right')   
    t.set_x(2.4)

typestr = ''
if use_scatter:
    typestr+='_zccode'
if only_failed:
    typestr+= '_only_failed'

if scatter_type == 'angle_location':
    #plt.ylim([0, 40])
    #plt.xlim([0, 10])
    print 'DEACTIVATED X/Y LIMS for TESTING!'

plt.xlabel(xlabel)
plt.ylabel(ylabel)

plt.suptitle(suptitle)
#plt.savefig('%s%s.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')

    
hax = fig.add_subplot(gs[1])
if only_failed:
    concat = num.concatenate((results_gotit, results_no))
else:
    concat = results
depths = set(concat.T[3])
print 'depth: ', depths
his = hax.hist(concat.T[3],
         len(depths), 
         orientation='horizontal',
         facecolor='gray', 
         normed=0,
         weights=num.zeros_like(concat.T[3])+1./concat.T[3].size,
         align='mid')
         
xmin, xmax = hax.get_xlim()
hax.set_xticks(num.arange(0,xmax%1,0.1))
hax.set_yticks(cticks)
hax.set_yticklabels([])
hax.set_ylim(vmin, vmax)
xticks = hax.get_xticklabels()
plt.setp(xticks, rotation=45)
hax.xaxis.grid(True, linestyle='-', which='major', color='grey', alpha=0.5)
fig.subplots_adjust(wspace=0.)

formatter = FuncFormatter(to_percent)

plt.gca().xaxis.set_major_formatter(formatter)
plt.xlabel('')

plt.savefig('%s%s.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')

if misfit is not None:
    figmisfit = plt.figure(figsize=(4,3), dpi=100) #, frameon=False, tight_layout=True)
    axmisfit= figmisfit.add_subplot(111)
    
    # Grid data:
    xg, yg = num.mgrid[X.min():X.max():100j, Y.min():Y.max():100j]

    if not nogrid:
        vg = griddata((X,Y), misfit, (xg,yg), method='cubic')

    sc = axmisfit.scatter(X,Y, c=misfit, s=8, lw=0.1, zorder=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar(sc, label='misfit M')
    plt.ylim([0, y_max])
    plt.xlim([0, x_max])
    plt.savefig('%s%s_misfit.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')
    #if not nogrid:
    #    plt.contourf(xg,yg,vg, zorder=0)


if scaling is not None:
    figscaling = plt.figure(figsize=(4,3), dpi=100) #, frameon=False, tight_layout=True)
    axscaling = figscaling.add_subplot(111)
    
    # Grid data:
    xg, yg = num.mgrid[X.min():X.max():100j, Y.min():Y.max():100j]

    if not nogrid:
        vg = griddata((X,Y), scaling, (xg,yg), method='cubic')

    sc = axscaling.scatter(X,Y, c=scaling, s=8, lw=0.1, zorder=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar(sc, label='scaling factor')
    plt.ylim([0, y_max])
    plt.xlim([0, x_max])
    plt.savefig('%s%s_scaling.pdf'%('.'.join(file_name.split('.')[:-1]), typestr), transparent=True, pad_inches=0.01, bbox_inches='tight')
    #if not nogrid:
    #    plt.contourf(xg,yg,vg, zorder=0)

plt.show()
