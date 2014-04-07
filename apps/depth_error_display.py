from guts import *
from derec.core import yamlTrace, TestCase, TestCaseData, TestCaseSetup
from derec.optics import OpticBase 
import matplotlib.pyplot as plt
from collections import defaultdict
from derec import derec_utils as du

fn = 'results/depth_error_4.0.yaml'
f = open(fn, 'r')
data1 = load_string(f.read())
f.close()

lines1 = TestCase.lines_dict(data1.candidates)
lines2 = TestCase.lines_dict(data1.references)
lines3 = TestCase.lines_dict(data1.processed_candidates)
lines4 = TestCase.lines_dict(data1.processed_references)

lines_resort = defaultdict(dict)
traces_resort = defaultdict(dict)

for s1,t1,l1 in TestCase.iter_dict(data1.candidates):
    for s2, t2, l2 in TestCase.iter_dict(data1.processed_candidates):
        if s1.__dict__==s2.__dict__ and t1.__dict__==t2.__dict__:
            traces_match = [du.yamlTrace2pyrockoTrace(l1), du.yamlTrace2pyrockoTrace(l2)]
            traces_match[0].set_codes('a','','','')
            traces_match[1].set_codes('b','','','')
            traces_resort[s1][t1] = traces_match 

for s1,t1,l1 in TestCase.iter_dict(lines1):
    for s2, t2, l2 in TestCase.iter_dict(lines2):
        for s3, t3, l3 in TestCase.iter_dict(lines3):
            for s4, t4, l4 in TestCase.iter_dict(lines4):

                if s1.__dict__==s2.__dict__==s3.__dict__==s4.__dict__ and t1.__dict__==t2.__dict__==t3.__dict__==t4.__dict__:
                    l3.set_linestyle('--')
                    l4.set_linestyle('--')
                    l1.set_color('r')
                    l3.set_color('r')
                    lines_match = [l1, l2, l3, l4]
                    lines_resort[s1][t1] = lines_match

fig = plt.figure()
ax = fig.add_subplot(111)
[ax.add_line(L) for L in lines_match]

ax.autoscale()
#plt.plot(lines2.values()[0].values()[0])

#opticbase = OpticBase(data)
#opticbase.waveforms_plot()
#opticbase.plot_marker_vs_distance()
#opticbase.plot_z_components(data.candidates, data.candidates_markers)


# Plot rise time vs. depth
plt.show()

