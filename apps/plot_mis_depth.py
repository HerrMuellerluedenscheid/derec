from guts import *
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
# change the *depth* to see a different depth.
#-------------------------------------------------------==


depth = 2000.
name = 'regional'
descriptor = ''

if True:
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
                extra=str(optic.test_case_setup.test_parameter_value)+\
                                '.'.join(k))
        v.savefig(fn)

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
    fn = file_name_path(name, 
            test_parameter=test_parameter,
            descriptor='', 
            test_type='z_diff_vs_tp.pdf', 
            extra='.'.join(k))
    fig.savefig(fn)
    #plt.show()


#............................................................................
# stack plot of one(!) result from different depths:
# reference traces as grey shaded areas
#..............................................
#colormap = plt.get_cmap('Blues')
#def scalez255(val):
#    return (val-min(depths))/max(depths)*255

if False:
    depths=[2000.,3000.]
    for k, opt in optics.iteritems():
        test_parameter = opt.test_case_setup.test_parameter 
        fig = plt.figure()
        ax = opt.stack_plot(depths=depths)
        #for a in ax.values():
        #    for l in a.get_lines():
        #        l.set_color(colormap(scalez255(float(\
        #            l.get_label().split()[0]))))

        fn = file_name_path(name, 
                test_parameter=test_parameter,
                descriptor='', 
                test_type='stack_z%s.pdf'%('_'.join(map(str, depths))), 
                extra=opt.test_case_setup.test_parameter_value)

        fig.savefig(fn)
        #plt.show()

#...............................................
# plot z components of chosen optic

#buggy
if False:
    opt.plot_z_components(traces_dict=data[0].processed_candidates,
            markers_dict=data[0].candidates_markers ,
            sources=opt.get_sources_where({'depth':depth}))
    plt.show()

