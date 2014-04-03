from pyrocko.gf import *
from pyrocko import model, gui_util, pile, trace, moment_tensor, io
from guts import *
from matplotlib import pyplot as plt

import progressbar
import os
from derec import derec_utils as du
from derec.core import *
import numpy as num


pjoin = os.path.join
km = 1000.

if __name__ ==  "__main__":

    fn = 'sample_test_case_setup.yaml'
    test_case_setup = load_string(open(fn,'r').read())

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]

    # load stations from file:
    stations = model.load_stations(pjoin(selfdir,
                            '../reference_stations_castor_selection.txt'))

    zoffset = 2000
    depths=num.linspace(test_case_setup.reference_source.depth-zoffset, 
                        test_case_setup.reference_source.depth+zoffset, 
                        9)

    print depths, '<- depths'

    reference_request = make_reference_trace(test_case_setup.reference_source,
                                             test_case_setup.targets, 
                                             test_case_setup.engine)

    reference_seismograms = du.response_to_dict(reference_request)
    test_case_dict = {}

    rise_times = num.linspace(1.0,10,4)
    for i, rise_time in enumerate(rise_times):
        print "%s of %s" % (i+1, len(rise_times))

        stf = [[0,rise_time], [0,1]]

        test_case = TestCase( test_case_setup )

        for tr in TestCase.iter_dict(reference_seismograms, only_values=True):
            du.add_random_noise_to_trace(tr, A=0.00001)

        test_case.set_raw_references(reference_seismograms)

        # considering that these traces are 'real' traces. Thus, 
        # stf needs to be applied to raw traces.
        test_case.raw_references = du.apply_stf(test_case.raw_references, 
                                test_case_setup.source_time_function)

        extended_ref_marker = du.chop_ranges(test_case_setup.reference_source, 
                                            test_case_setup.targets, 
                                            test_case.store,
                                            test_case_setup.phase_ids_start,
                                            perc=1.0,
                                            t_shift_frac=0.3)

        test_case.set_reference_markers(extended_ref_marker)

        D = Doer(test_case)
        test_case_dict[rise_time] = test_case
        test_case.yaml_dump(fn='results/depth_error_%s.yaml'%rise_time)

    fig = plt.figure()


    for key, test_case in test_case_dict.items():
        
        true_z = test_case.test_case_setup.reference_source.depth
        best_source, best_fit = test_case.get_lowest_misfit_data()
        zdiff = true_z - best_source.depth
        print key, zdiff
        plt.plot(key, zdiff, 'o')

    plt.show()
 


