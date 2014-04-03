from guts import *
from derec.core import yamlTrace, TestCase, TestCaseData, TestCaseSetup
from derec.optics import OpticBase 
import matplotlib.pyplot as plt

fn = 'results/depth_error_4.0.yaml'
f = open(fn, 'r')
data = load_string(f.read())
f.close()
opticbase = OpticBase(data)
opticbase.waveforms_plot()
opticbase.plot_marker_vs_distance()
opticbase.plot_z_components(data.candidates, data.candidates_markers)
plt.show()

