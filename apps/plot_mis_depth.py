from guts import *
from derec.core import TestCase, TestCaseData, TestCaseSetup
from derec.optics import OpticBase 
import matplotlib.pyplot as plt
from collections import defaultdict
from derec import derec_utils as du
import progressbar
import glob
import sys

def add_line_to_ca(line):
    line.set_figure(plt.gcf())
    plt.gca().add_line(line)

def combinator(lists_dict):
    """Merge dicts into dicts with lists. """
    result = defaultdict(list)
    for d in lists_dict:
        for k,v in d.iteritems():
            result[k].append(v)
    return result

depth = 2000.

file_names = sys.argv[1:]
if not file_names:
    print '''usage: python plot_mis_depth.py 'files_of_interest' '''
    sys.exit(0)

pbar = progressbar.ProgressBar(len(file_names))

optics = {}
data = []
for i, fn in enumerate(file_names):
    pbar.update(i)
    #f = open(fn, 'r')
    for date in load_all(filename=fn):
        #f.close()
        data.append(date)
        ob = OpticBase(date)
        optics.update({ob.test_case_setup.test_parameter_value: ob})

pbar.finish()

sorted_keys = sorted(optics.keys())
all_targets = optics.values()[0].targets

#.........................................................................
# Plot: 1 per station, subplots: one per optics.
# show reference trace plotted over best fitting candidate
#-------------------------------------------------------==
if False:
    fig_dict = OpticBase.figure_dict([t.codes for t in all_targets]) 
    for i,k in enumerate(sorted_keys):
        optic = optics[k]
        for source in optic.get_sources_where({'depth':depth}):
            for target in optic.distance_sort_targets(source):
                fig = plt.figure(fig_dict[target.codes].number)
                ax = fig.add_subplot(len(optics), 1, i)
                l1 = optic.get_processed_reference_line(source, target)
                l2 = optic.get_processed_candidate_line(source, target)
                map(add_line_to_ca, [l1,l2])
                ax.autoscale()
    plt.show()

#..........................................................................
# Plot depth difference versus test_parameter:
#---------------------------------------------
if False:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for date in data:
        setup = date.test_case_setup
        test_parameter = setup.test_parameter 
        test_parameter_value = setup.test_parameter_value
        best_s = min(date.misfits, key=date.misfits.get)
        zdiff = date.test_case_setup.reference_source.depth-best_s.depth
        ax.plot(test_parameter_value, zdiff, 'o')
    ax.autoscale()
    ax.set_xlabel(test_parameter)
    plt.show()

#............................................................................
# stack plot:
# looks pretty much the same as first plot. 
#..............................................
opt = optics.values()[0]
if True:
    ax1 = opt.stack_plot( markers_dict=data[0].candidates_markers, depths=[depth])

#...............................................
# plot z components of chosen optic

if True:
    opt.plot_z_components(traces_dict=data[0].processed_candidates,
            markers_dict=data[0].candidates_markers ,
            sources=opt.get_sources_where({'depth':depth}))
    plt.show()

