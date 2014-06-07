import copy
import progressbar
import os
import numpy as num
import copy
from matplotlib import pyplot as plt

from pyrocko.gf import *
from pyrocko import model, gui_util, pile, trace, moment_tensor, io, parimap
from pyrocko.guts import *
from derec.optics import OpticBase
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
        print i, 'of', len(test_parameter_values)

        setattr(test_case_setup, 'test_parameter', test_parameter)

        test_case_setup.test_parameter_value = float(parameter_value) 

        reset_source(reference_source, __reference_source_copy)

        if test_case_setup.test_parameter=='source_time_function':
            stf[0][1] = float(parameter_value)
            setattr(test_case_setup, test_parameter, stf)
        else:
            setattr(test_case_setup, 'source_time_function', __stf)
            setattr(reference_source, test_parameter, float(parameter_value))

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
        print 'start'
        test_case.process()
        print 'done'
        if debug:
            op = OpticBase(test_case)
            op.stack_plot()
            plt.show() 
            import pdb
            pdb.set_trace()

        
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
    #name = 'rotenburg_standard' 
    #choice = 'rotenburg'
    name = 'doctar_fewerdepths'
    choice = 'doctar'
    #name = 'regional_bandpass' 
    #name = 'global'
    #name = 'castor'
    #choice = 'castor'

    add_noise = False
    test_model = False
    debug = False

    descriptor = ''
    description = '''noise free test. This time with bandpass. longer static length 
    to see if source time function dependence gets better'''


    if choice=='doctar':

        norm = 2 
        taper = trace.CosFader(xfade=1.)  
         
        z, p, k = butter(2, [0.7*num.pi*2, 6.0*num.pi*2.],  
                           'band',  
                           analog=True,  
                           output='zpk') 
         
        z = map(complex, z) 
        p = map(complex, p) 
        k = complex(k) 
         
        fresponse = trace.PoleZeroResponse(z,p,k) 
        fresponse.validate() 
     
        misfit_setup = trace.MisfitSetup(norm=norm, 
                                         taper=taper, 
                                         domain='time_domain', 
                                         filter=fresponse) 


        fn = 'test_case_setup.yaml'
        #ok:
        store_id = 'doctar_mainland_20Hz'
        #test_case_setup = load_string(open(fn,'r').read())
        store_dirs = [pjoin(derec_home, 'fomostos'), '/scratch/local1/marius']

        engine = LocalEngine(store_superdirs=store_dirs)    

        reference_event = model.Event(load=pjoin(derec_home, \
                'mseeds/doctar/doctar_2011-11-01/doctar_2011-11-01_quakefile.dat'))
        reference_source = DCSource.from_pyrocko_event(\
                reference_event)
        stations = model.load_stations(pjoin(derec_home, 'mseeds/doctar/doctar_2011-11-01'\
                , 'stations.txt'))

        targets = []
        for s in stations:
            targets.extend(Target.from_pyrocko_station(s, store_id=store_id))
        

        #num_depths = 11
        #offset = 2000
        #depths = num.linspace(reference_source.depth-offset,
        #        reference_source.depth+offset, num_depths)
        depths = num.array(range(1400, 10000, 200))
        print 'ref depth', reference_source.depth
        print depths
        test_sources = du.test_event_generator(reference_source, depths)
        
        stf = [[0., 0.2],[0.,1.]]
            
        phase_ids_start = ['p','P']

        test_case_setup = TestCaseSetup(reference_source=reference_source,
                                       sources=test_sources,
                                       targets=targets,
                                       engine=engine, 
                                       store_id=store_id,
                                       misfit_setup=misfit_setup,
                                       source_time_function=stf,
                                       number_of_time_shifts=51,
                                       percentage_of_shift=30.,
                                       phase_ids_start=phase_ids_start,
                                       static_length=2.8,
                                       marker_perc_length=5.0,
                                       marker_shift_frac=0.3,
                                       depths=depths) 

    if choice=='rotenburg':
        fn = 'test_case_setup.yaml'
        test_case_setup = load_string(open(fn,'r').read())
        store_id = 'crust2_dd'
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
        norm = 2 
        taper = trace.CosFader(xfade=1.)  
         
        z, p, k = butter(2, [0.7*num.pi*2, 6.0*num.pi*2.],  
                           'band',  
                           analog=True,  
                           output='zpk') 
         
        z = map(complex, z) 
        p = map(complex, p) 
        k = complex(k) 
         
        fresponse = trace.PoleZeroResponse(z,p,k) 
        fresponse.validate() 
     
        misfit_setup = trace.MisfitSetup(norm=norm, 
                                         taper=taper, 
                                         domain='time_domain', 
                                         filter=fresponse) 

        #fn = 'test_case_setup_castor.yaml'
        test_case_setup = load_string(open(fn,'r').read())
        store_id = test_case_setup.store_id
        reference_source = test_case_setup.reference_source
        engine = test_case_setup.engine

        store_id = 'castor_20Hz'
        #test_case_setup = load_string(open(fn,'r').read())
        store_dirs = [pjoin(derec_home, 'fomostos'), '/scratch/local1/marius']

        engine = LocalEngine(store_superdirs=store_dirs)    

        reference_event = model.Event(load=pjoin(derec_home, \
                'mseeds/castor/castor_event_2013-10-01.dat'))
        reference_source = DCSource.from_pyrocko_event(\
                reference_event)
        stations = model.load_stations(pjoin(derec_home, 'mseeds/castor'\
                , 'stations.txt'))

        targets = []
        for s in stations:
            targets.extend(Target.from_pyrocko_station(s, store_id=store_id))
        
        num_depths = 11
        offset = 2000
        depths = num.linspace(reference_source.depth-offset,
                reference_source.depth+offset, num_depths)
        test_sources = du.test_event_generator(reference_source, depths)
        
        stf = [[0., 1.],[0.,1.]]
            
        phase_ids_start = ['p','P']
        test_case_setup = TestCaseSetup(reference_source=reference_source,
                                       sources=test_sources,
                                       targets=targets,
                                       engine=engine, 
                                       store_id=store_id,
                                       misfit_setup=misfit_setup,
                                       source_time_function=stf,
                                       number_of_time_shifts=21,
                                       percentage_of_shift=30.,
                                       phase_ids_start=phase_ids_start,
                                       static_length=3.5,
                                       marker_perc_length=5.0,
                                       marker_shift_frac=0.55,
                                       depths=depths) 

    ref_lat = reference_source.lat
    ref_lon = reference_source.lon 
    ref_depth = reference_source.depth
    ref_strike = reference_source.strike
    ref_dip = reference_source.dip
    ref_rake = reference_source.rake
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
    #if robust_check:
    #    targets = test_case_setup.targets
    #    target_on_source = Target(lat=ref_lat,
    #                              lon=ref_lon,
    #                              codes=('','','','BHZ'))

    #    num_stations = len()
    #    test_case_setup.targets = setup_targets(target_on_source,
    #                                            num_stations,
    #                                            100*km)

    #test_case_setup.number_of_time_shifts = 21
    #test_case_setup.static_length = 5.5

    __reference_source_copy = du.clone(reference_source)

    zoffset = 2000

    #depths=num.linspace(reference_source.depth-zoffset, 
    #                    reference_source.depth+zoffset, 
    #                    11)


    test_case_setup.engine.store_superdirs = store_dirs
    stf = [[0., 0.2], [0.,1.]]
    __stf = copy.deepcopy(stf)

    noise = []
    if add_noise:
        for tr in io.load(glob.glob(pjoin(derec_home, 'mseeds', 'iris_data',
            'checked_noise')+'/*')):
            noise.extend(tr)
    
    reference_seismograms = make_reference_trace(reference_source,
                                             test_case_setup.targets, 
                                             test_case_setup.engine,
                                             stf,
                                             noise)

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

    # note: doctar: 0.1-0.5 balken bei richtier tiefe...
    test_parameter_values = [num.linspace(0.02, 1.2, 31),
                             num.linspace(ref_strike-45., ref_strike+45., 15),
                             num.linspace(ref_dip-45., ref_dip+45., 15),
                             num.linspace(ref_rake-45., ref_rake+45., 15),
                             num.linspace(lat_shift_n, lat_shift_p, 15),
                             num.linspace(lon_shift_n, lon_shift_p, 15)]

    for i, tpset in enumerate(zip(test_parameter, test_parameter_values)):
        print i+1, 'of', len(test_parameter)
        do_run(tpset) 
    print 'finised %s'% name
