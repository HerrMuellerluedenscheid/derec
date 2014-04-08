from guts import *
from derec.core import yamlTrace, TestCase, TestCaseData, TestCaseSetup
from derec.optics import OpticBase 
import matplotlib.pyplot as plt
from collections import defaultdict
from derec import derec_utils as du
import glob

file_names = glob.glob('results/*yaml')

data = []
for fn in file_names:
    f = open(fn, 'r')
    data.append(load_string(f.read()))
    f.close()

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

