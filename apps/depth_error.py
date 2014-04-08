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

    zoffset = 1000

    depths=num.linspace(test_case_setup.reference_source.depth-zoffset, 
                        test_case_setup.reference_source.depth+zoffset, 
                        3)

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

    reference_seismograms = du.response_to_dict(reference_request)
    test_case_dict = {}

    test_case = TestCase( test_case_setup )

    reference_seismograms = du.apply_stf(reference_seismograms, 
                            test_case_setup.source_time_function)

    test_case.set_raw_references(reference_seismograms)

    extended_ref_marker = du.chop_ranges(test_case.reference_source, 
                                test_case.targets, 
                                test_case.store,
                                test_case.phase_ids_start,
                                perc=1.0,
                                t_shift_frac=0.3,
                                use_cake=True)

    test_case.set_reference_markers(extended_ref_marker)

    D = Doer(test_case)
    import pdb
    pdb.set_trace()
    print test_case.misfits
    test_case.yaml_dump(fn='results/depth_error_1.yaml')

