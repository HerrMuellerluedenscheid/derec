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

    fn = 'sample_test_case_setup.yaml'
    test_case_setup = load_string(open(fn,'r').read())

    zoffset = 2000

    depths=num.linspace(test_case_setup.reference_source.depth-zoffset, 
                        test_case_setup.reference_source.depth+zoffset, 
                        21)

    # overwriting sources:
    test_case_setup.sources = [DCSource(lat=test_case_setup.reference_source.lat,
                            lon=test_case_setup.reference_source.lon,
                            depth=float(depth),
                            time=test_case_setup.reference_source.time,
                            strike=test_case_setup.reference_source.strike,
                            dip=test_case_setup.reference_source.dip,
                            rake=test_case_setup.reference_source.rake,
                            magnitude=test_case_setup.reference_source.magnitude) 
                                          for depth in depths]

    print depths, '<- depths'
    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [derec_home + '/fomostos']
    test_case_setup.engine.store_superdirs = store_dirs
    reference_request = make_reference_trace(test_case_setup.reference_source,
                                             test_case_setup.targets, 
                                             test_case_setup.engine)


    rise_times = num.linspace(0.5,4.5,30)
    for rise_time in rise_times:
        test_case = TestCase( test_case_setup )
        test_case_setup.test_parameter = 'rise_time'
        test_case_setup.test_parameter_value = float(rise_time)

        # overwriting reference sources' stf.
        # the stf of the candidates is taken from the setup and therefore stays
        # the same.
        stf = [[0., float(rise_time)], [0.,1.]]
        reference_seismograms = du.response_to_dict(reference_request)
        test_case_dict = {}
        reference_seismograms = du.apply_stf(reference_seismograms, 
                                                        stf)

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
        test_case.yaml_dump(fn='results/depth_error_%s.yaml'%rise_time)
        #from depth_error_display import make_compare_plots as mcp
        #mcp(test_case.candidates, test_case.references,
        #        test_case.processed_candidates, test_case.processed_references)


