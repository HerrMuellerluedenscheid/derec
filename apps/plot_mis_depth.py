from guts import *
from derec.core import TestCase, TestCaseData, TestCaseSetup
from derec.optics import OpticBase 
import matplotlib.pyplot as plt
from collections import defaultdict
from derec import derec_utils as du
import progressbar
import glob
import sys

file_names = sys.argv[1:]

pbar = progressbar.ProgressBar(len(file_names))

data = []
for i, fn in enumerate(file_names):
    pbar.update(i)
    f = open(fn, 'r')
    data.append(load_string(f.read()))
    f.close()

pbar.finish()
#..........................................................................
# Plot depth difference versus test_parameter:
#---------------------------------------------
#fig = plt.figure()
#ax = fig.add_subplot(111)
#for date in data:
#    setup = date.test_case_setup
#    test_parameter = setup.test_parameter 
#    test_parameter_value = setup.test_parameter_value
#    best_s = min(date.misfits, key=date.misfits.get)
#    zdiff = date.test_case_setup.reference_source.depth-best_s.depth
#    ax.plot(test_parameter_value, zdiff, 'o')
#ax.autoscale()
#ax.set_xlabel(test_parameter)
#plt.show()
#............................................................................

optics = OpticBase(data[0])
optics.plot_z_components(traces_dict=data[0].processed_candidates,
        markers_dict=data[0].candidates_markers ,
        sources=optics.get_sources_where({'depth':2000}))
plt.show()
