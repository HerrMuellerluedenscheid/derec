from pyrocko.guts import *
from derec.core import yamlTrace, TestCase, TestCaseData, TestCaseSetup
from derec.optics import OpticBase 
import matplotlib.pyplot as plt
from collections import defaultdict
from derec import derec_utils as du
import glob

#fn = 'results/depth_error_1.yaml'
#f = open(fn, 'r')
#data1 = load_string(f.read())
#f.close()
#
#ttt=data1.candidates.values()[0].values()[0]
#print ttt.get_xdata()
#print type(ttt)
#
#lines1 = TestCase.lines_dict(data1.candidates)
#lines2 = TestCase.lines_dict(data1.references)
#lines3 = TestCase.lines_dict(data1.processed_candidates)
#lines4 = TestCase.lines_dict(data1.processed_references)
#

#for s1,t1,l1 in TestCase.iter_dict(lines1):
#    for s2, t2, l2 in TestCase.iter_dict(lines2):
#        for s3, t3, l3 in TestCase.iter_dict(lines3):
#            for s4, t4, l4 in TestCase.iter_dict(lines4):

def make_compare_plots(d1, d2, d3, d4):
    for s1,t1,l1 in TestCase.iter_dict(d1):
        for s2, t2, l2 in TestCase.iter_dict(d2):
            for s3, t3, l3 in TestCase.iter_dict(d3):
                for s4, t4, l4 in TestCase.iter_dict(d4):

                    if s1.__dict__==s2.__dict__==s3.__dict__==s4.__dict__ and t1.__dict__==t2.__dict__==t3.__dict__==t4.__dict__:


                        lines_match = [l1, l2, l3, l4]
                        ptrac = [du.yamlTrace2pyrockoTrace(l) for l in lines_match]
                        [ti.set_codes(str(i),'','','') for i,ti in enumerate(ptrac)]
                        
                        fig = plt.figure()
                        ax = fig.add_subplot(111)
                        [ax.plot(l.get_xdata(), l.get_ydata()) for l in ptrac]
                        ax.autoscale()



    #plt.plot(lines2.values()[0].values()[0])
    #plt.figure()

    #opticbase = OpticBase(data1)
    #opticbase.waveforms_plot()
    #opticbase.plot_marker_vs_distance()
    #opticbase.plot_z_components(data1.candidates, data1.candidates_markers)
    #opticbase.process_compare(ignore=['depth'])


    # Plot rise time vs. depth
    plt.show()

