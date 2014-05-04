import matplotlib
matplotlib.use('Agg')
from pyrocko.guts import *
from derec.core import TestCase, TestCaseData, TestCaseSetup
from derec.optics import OpticBase, gca_label 

import matplotlib.pyplot as plt
from collections import defaultdict
from derec import derec_utils as du
import progressbar
import glob
import sys
import os

pjoin = os.path.join

def add_line_to_ca(line):
    line.set_figure(plt.gcf())
    plt.gca().add_line(line)

def file_name_path(name, test_parameter, descriptor, test_type, extra ):
    fn = '%s/%s%s/%s_%s.pdf'%(name,
                             test_parameter,
                             descriptor,
                             test_type,
                             extra) 
    return fn

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
    name = fn.split('/')[0]
    for date in load_all(filename=fn):
        data.append(date)
        ob = OpticBase(test_case_data=date)
        optics.update({ob.test_case_setup.test_parameter_value: ob})

pbar.finish()

sorted_keys = sorted(optics.keys())
all_targets = optics.values()[0].targets

#.........................................................................
# Plot: 1 per station, subplots: one per optics.
# compares the influence of different parameters on waveforms
# show reference trace plotted over best fitting candidate
# change the *depth* to see a different depth.
#-------------------------------------------------------==

depth = 2000.
#name = 'regional_new'
descriptor = ''

if False:
    print '(1) start generating processed plots'
    fig_dict = OpticBase.figure_dict([t.codes for t in all_targets]) 
    for i,k in enumerate(sorted_keys):
        optic = optics[k]
        for source in optic.get_sources_where({'depth':depth}):
            for target in optic.distance_sort_targets(source):
                fig = plt.figure(fig_dict[target.codes].number)

                ax = fig.add_subplot(len(optics), 1, i)
                l1 = optic.get_processed_reference_line(source, target)
                l2 = optic.get_processed_candidate_line(source, target)
                l2.set_color((1,0,0,1))
                l2.set_label('cand')
                l1.set_label('ref')
                map(add_line_to_ca, [l1,l2])
                ax.autoscale()
                gca_label(optic.test_case_setup.test_parameter_value)
                plt.suptitle('.'.join(target.codes))
        plt.legend()

    for k, v in fig_dict.iteritems():
        v.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fn=file_name_path(name, 
                test_parameter=optic.test_case_setup.test_parameter,
                descriptor='', 
                test_type='wvfrm_ref_vs_best_2000m', 
                extra='')
        v.savefig(fn)
    print 'done.'

#..........................................................................
# Plot depth difference versus test_parameter:
#---------------------------------------------
if False:
    print '(2)start generating diff plot'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for date in data:
        setup = date.test_case_setup
        test_parameter = setup.test_parameter 
        test_parameter_value = setup.test_parameter_value
        best_s = min(date.misfits, key=date.misfits.get)
        test_parameter_value = setup.test_parameter_value
        best_s = min(date.misfits, key=date.misfits.get)
        zdiff = date.test_case_setup.reference_source.depth-best_s.depth
        ax.plot(test_parameter_value, zdiff, 'o')
    ax.autoscale()
    ax.set_xlabel(test_parameter)
    fn = file_name_path(name, 
            test_parameter=test_parameter,
            descriptor='', 
            test_type='z_diff_vs_tp.pdf', 
            extra='')
    fig.savefig(fn)
    #plt.show()
    print 'done'


#............................................................................
# stack plot of one(!) result from different depths:
# reference traces as grey shaded areas
#..............................................

if True:
    print '(3)start generating stack plot'
    depths=[1000.,2000.,3000.]
    for k, opt in optics.iteritems():
        fig = plt.figure()
        axs = opt.stack_plot(depths=depths)#, fig=fig)
        test_parameter = opt.test_case_setup.test_parameter 
        fn = file_name_path(name, 
                test_parameter=test_parameter,
                descriptor='', 
                test_type='stack', 
                extra=opt.test_case_setup.test_parameter_value)

        fig.savefig(fn)
    print 'done'
#...............................................
# plot z components of chosen optic

if False:
    channel = 'Z'
    for param_value, opt in optics.items():
        test_parameter = opt.test_case_setup.test_parameter 
        srcs = opt.get_sources_where({'depth':depth})
        for src, fig in opt.plot_channel(traces_dict=opt.candidates,
                                markers_dict=opt.test_case_data.candidates_markers,
                                channel=channel,
                                sources=srcs):

            fn = file_name_path(name, 
                    test_parameter=test_parameter,
                    descriptor='', 
                    test_type='channel_%s_%s'%(src.depth, channel),
                    extra=param_value)

            fig.savefig(fn)

#-----------------------------------------------
# Plot the misfit dict
if False:
    for param_value, opt in optics.items():
        fig = opt.plot_misfits() 
        fn = file_name_path(name, 
                test_parameter=opt.test_case_setup.test_parameter,
                descriptor='', 
                test_type='misfits_',
                extra=param_value)
        fig.savefig(fn)
