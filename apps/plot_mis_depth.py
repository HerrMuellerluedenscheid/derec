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
    rise_time = date.test_case_setup.source_time_function
    z = min(date.misfits, key=date.misfits.get)
    zdiff = date.test_case_setup.reference_source.depth-z
    print z, rise_time
    ax.plot(rise_time, zdiff)

# Plot rise time vs. depth
plt.show()

