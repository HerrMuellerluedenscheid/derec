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


def setup_targets(ref_target, num_stations, field_range, chas=['BHE','BHN','BHZ']):
    """
    Setup randomized targets.
    """
    e_shifts = num.random.uniform(-field_range, field_range, num_stations)
    n_shifts = num.random.uniform(-field_range, field_range, num_stations)
    return [Target(lat=ref_target.lat,
                   lon=ref_target.lon,
                   depth=ref_target.depth,
                   codes=(ref_target.codes[0],
                          '%s_'%ti+ref_target.codes[1],
                          '',
                          c),
                   north_shift = n_shifts[ti],
                   east_shift = e_shift)
                      for ti,e_shift in enumerate(e_shifts) 
                      for c in chas]




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
        test_case.process()
        
        base_dir = pjoin(name, test_case_setup.test_parameter, descriptor)
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)

        test_case.yaml_dump(fn='%s/%s%s/depth_error_%s.yaml'%(name,
                                    test_case_setup.test_parameter,
                                    descriptor,
                                    test_case_setup.test_parameter_value))

if __name__ ==  "__main__":
    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [derec_home + '/fomostos']

    #name = 'local'
    name = 'rotenburg' 
    choice='rotenburg'
    #name = 'regional_bandpass' 
    #name = 'global'

    descriptor = ''
    description = 'noise free test. This time with bandpass.'

    fn = 'test_case_setup.yaml'
    test_case_setup = load_string(open(fn,'r').read())
    store_id = test_case_setup.store_id

    if choice=='rotenburg':
        reference_event = model.Event(load=pjoin(derec_home, \
                'mseeds/RotenburgIris/rotenburg_quakefile.dat'))
        reference_source = DCSource.from_pyrocko_event(\
                reference_event)
        test_case_setup.reference_source = reference_source
        stations = model.load_stations(pjoin(derec_home, 'mseeds/RotenburgIris'\
                , 'stations.txt'))
        targets = []
        for s in stations:
            targets.extend(Target.from_pyrocko_station(s, store_id=store_id))
        test_case_setup.targets = targets

    elif choice=='xizang':
        pass
        # filter [0.05-1.5 Hz]
    elif choice=='castor':
        reference_source = test_case_setup.reference_source

    ref_lat = reference_source.lat
    ref_lon = reference_source.lon 
    ref_depth = reference_source.depth
    ref_strike = reference_source.strike
    ref_dip = reference_source.dip
    ref_rake = reference_source.rake
    ref_id = reference_source.store_id
    #z, p, k = butter(2, [0.001*num.pi*2, 1.0*num.pi*2.],  
    #                   'band',  
    #                   analog=True,  
    #                   output='zpk') 
    # 
    #z = map(complex, z) 
    #p = map(complex, p) 
    #k = complex(k) 
    #fresponse = trace.PoleZeroResponse(z,p,k)
    #fresponse.validate()
    #test_case_setup.misfit_setup.filter = fresponse
    if robust_check:
        targets = test_case_setup.targets
        target_on_source = Target(lat=ref_lat,
                                  lon=ref_lon,
                                  codes=('','','','BHZ'))

        num_stations = len()
        test_case_setup.targets = setup_targets(target_on_source,
                                                num_stations,
                                                100*km)

    test_case_setup.number_of_time_shifts = 21
    test_case_setup.static_length = 5.

    __reference_source_copy = reference_source.clone()

    zoffset = 1000

    depths=num.linspace(reference_source.depth-zoffset, 
                        reference_source.depth, 
                        21)

    depths = [float(d) for d in depths]
    print depths, '<- depths'
    test_case_setup.depths = depths

    test_case_setup.engine.store_superdirs = store_dirs
    stf = [[0., 1.], [0.,1.]]
    __stf = copy.deepcopy(stf)

    reference_seismograms = make_reference_trace(reference_source,
                                             test_case_setup.targets, 
                                             test_case_setup.engine,
                                             stf, 
                                             noisedir = pjoin(derec_home,
                                                              'mseeds',
                                                              'iris_data',
                                                              'restitute'))
    test_parameter = ['source_time_function',
                      'strike',
                      'dip',
                      'rake',
                      'latitude',
                      'longitude' ]
    if test_model:
        test_parameter+='store_id'


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

    test_parameter_values = [num.linspace(0.5, 5., 21),
                             num.linspace(ref_strike-45., ref_strike+45., 15),
                             num.linspace(ref_dip-45., ref_dip+45., 15),
                             num.linspace(ref_rake-45., ref_rake+45., 15),
                             num.linspace(lat_shift_n, lat_shift_p, 15),
                             num.linspace(lon_shift_n, lon_shift_p, 15)]

    for i, tpset in enumerate(zip(test_parameter, test_parameter_values)):
        print i+1, 'of', len(test_parameter)
        do_run(tpset) 


    
