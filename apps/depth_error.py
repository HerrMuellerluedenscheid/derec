from pyrocko.gf import *
from pyrocko import model, gui_util, pile, trace, moment_tensor, io
from pyrocko.guts import *
from matplotlib import pyplot as plt

import progressbar
import os
from derec import derec_utils as du
from derec.core import *
import numpy as num


pjoin = os.path.join
km = 1000.

if __name__ ==  "__main__":

    name = 'regional'
    description = 'noise free test. This time using more shifts to see if the'+\
    'errors decrease.'

    fn = 'test_case_setup.yaml'
    test_case_setup = load_string(open(fn,'r').read())

    zoffset = 2000

    depths=num.linspace(test_case_setup.reference_source.depth-zoffset, 
                        test_case_setup.reference_source.depth+zoffset, 
                        41)

    depths = [float(d) for d in depths]
    print depths, '<- depths'
    test_case_setup.depths = depths

    
    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [derec_home + '/fomostos']
    test_case_setup.engine.store_superdirs = store_dirs
    reference_request = make_reference_trace(test_case_setup.reference_source,
                                             test_case_setup.targets, 
                                             test_case_setup.engine)

    stf = [[0., 1.0], [0.,1.]]
    reference_seismograms = du.response_to_dict(reference_request)
    reference_seismograms = du.apply_stf(reference_seismograms, 
                                                    stf)

    strikes = num.arange(100.,160., 2.)

    for i, strike in enumerate(strikes):
        test_case_setup.reference_source.strike = float(strike)

        # overwriting sources:
        test_case_setup.sources = du.test_event_generator(
                                test_case_setup.reference_source, depths)

        print '%s of %s'%(i+1, len(strikes))
        test_case = TestCase( test_case_setup )

        test_case_setup.test_parameter = 'strike'
        test_case_setup.test_parameter_value = float(strike)

        test_case_dict = {}
        test_case.set_raw_references(reference_seismograms)

        #perc ist die Ausdehnung des zeitmarkers. perc=1 heisst (tmin-t0)*perc,
        # also der Ausschnitt ist so lang wie die Zeit zwischen t0 und tmin.
        extended_ref_marker = du.chop_ranges(test_case.reference_source, 
                                    test_case.targets, 
                                    test_case.store,
                                    test_case.phase_ids_start,
                                    perc=1.0,
                                    t_shift_frac=0.3,
                                    use_cake=True)

        test_case.set_reference_markers(extended_ref_marker)

        D = Doer(test_case)
        #from depth_error_display import make_compare_plots as mcp
        #mcp(test_case.candidates, test_case.references,
        #        test_case.processed_candidates, test_case.processed_references)
        #plt.show()
        test_case.yaml_dump(fn='%s/%s%s/depth_error_%s.yaml'%(name,
            test_case_setup.test_parameter,'',
            test_case_setup.test_parameter_value))
