import copy
import progressbar
import os
import numpy as num
import copy
from matplotlib import pyplot as plt

from pyrocko.gf import *
from pyrocko import model, gui_util, pile, trace, moment_tensor, io, parimap
from pyrocko.guts import *
from derec import derec_utils as du
from derec.core import *

def reset_source(source, ref_source):
    source.__dict__ = ref_source.__dict__.copy()


pjoin = os.path.join
km = 1000.

def do_run(tpvalues):    
    test_parameter = tpvalues[0]
    test_parameter_values = tpvalues[1]
    for i, parameter_value in enumerate(test_parameter_values):

        setattr(test_case_setup, 'test_parameter', test_parameter)

        test_case_setup.test_parameter_value = float(parameter_value) 

        reset_source(reference_source, __reference_source_copy)

        if test_case_setup.test_parameter=='source_time_function':
            stf[0][1] = float(parameter_value)
            setattr(test_case_setup, test_parameter, stf)
        else:
            setattr(test_case_setup, 'source_time_function', __stf)
            setattr(reference_source, test_parameter, float(parameter_value))
            print __stf, 'check! should be [01][01]'

        # overwriting sources:
        test_case_setup.sources = du.test_event_generator(
                                                reference_source, depths)

        test_case = TestCase( test_case_setup )
        test_case.set_raw_references(reference_seismograms)
        test_case_dict = {}

        extended_ref_marker = du.chop_ranges(test_case.reference_source, 
                                    test_case.targets, 
                                    test_case.store,
                                    test_case.phase_ids_start,
                                    perc=test_case_setup.marker_perc_length,
                                    static_length=test_case_setup.static_length,
                                    t_shift_frac=test_case_setup.marker_shift_frac,
                                    use_cake=True)

        test_case.set_reference_markers(extended_ref_marker)
        D = Doer(test_case)
        
        base_dir = pjoin(name, test_case_setup.test_parameter, descriptor)
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)

        test_case.yaml_dump(fn='%s/%s%s/depth_error_%s.yaml'%(name,
                                    test_case_setup.test_parameter,
                                    descriptor,
                                    test_case_setup.test_parameter_value))


if __name__ ==  "__main__":

    #name = 'local'
    name = 'regional_test' 
    #name = 'global'

    descriptor = ''
    description = 'noise free test. This time using more shifts to see if the'+\
    'errors decrease.'

    fn = 'test_case_setup.yaml'
    test_case_setup = load_string(open(fn,'r').read())
    reference_source = test_case_setup.reference_source
    __reference_source_copy = copy.deepcopy(reference_source)

    zoffset = 2000

    depths=num.linspace(reference_source.depth-zoffset, 
                        reference_source.depth+zoffset, 
                        21)

    depths = [float(d) for d in depths]
    print depths, '<- depths'
    test_case_setup.depths = depths

    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [derec_home + '/fomostos']
    test_case_setup.engine.store_superdirs = store_dirs
    stf = [[0., 1.], [0.,1.]]
    __stf = copy.deepcopy(stf)

    reference_seismograms = make_reference_trace(reference_source,
                                             test_case_setup.targets, 
                                             test_case_setup.engine,
                                             stf)

    test_parameter = ['source_time_function',
                      'strike',
                      'dip',
                      'rake',
                      'latitude',
                      'longitude' ]

    ref_lat = reference_source.lat
    ref_lon = reference_source.lon 
    ref_depth = reference_source.depth
    ref_strike= reference_source.strike
    ref_dip =  reference_source.dip
    ref_rake= reference_source.rake

    n_shift = 2000.
    e_shift = 2000.
    lon_shift_p = du.lat_lon_relative_shift(ref_lat, 
                                          ref_lon,
                                          east_shift=e_shift)

    lat_shift_p = du.lat_lon_relative_shift(ref_lat, 
                                          ref_lon,
                                          north_shift=n_shift)

    lon_shift_n = du.lat_lon_relative_shift(ref_lat, 
                                          ref_lon,
                                          east_shift=-e_shift)

    lat_shift_n = du.lat_lon_relative_shift(ref_lat, 
                                          ref_lon,
                                          north_shift=-n_shift)

    test_parameter_values = [num.linspace(0.5, 5., 15),
                             num.linspace(ref_strike-45., ref_strike+45., 15),
                             num.linspace(ref_dip-45., ref_dip+45., 15),
                             num.linspace(ref_rake-45., ref_rake+45., 15),
                             num.linspace(lat_shift_n, lat_shift_p, 15),
                             num.linspace(lon_shift_n, lon_shift_p, 15)]

    for i, tpset in enumerate(zip(test_parameter, test_parameter_values)):
        print i+1, 'of', len(test_parameter)
        do_run(tpset) 
    
