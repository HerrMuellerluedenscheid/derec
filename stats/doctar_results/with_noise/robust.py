from collections import defaultdict
import matplotlib.lines as pltlines
import matplotlib.pyplot as plt
import progressbar
import os
from derec import derec_utils as du
from derec import core
import numpy as num
import glob
from derec.yaml_derec import *
from derec import optics
from pyrocko.gf import *
from pyrocko import model, gui_util, trace, moment_tensor, io
from scipy.signal import butter
from scipy.ndimage import zoom
from pyrocko.guts import Object, Float, Int, String, Complex, Tuple, List, load_string, Dict
from pyrocko.guts_array import Array
import time
import progressbar


pjoin = os.path.join
km = 1000.

def pbar(i, num_tests, pb=None ):
    try:
        pb.update(i)
    except AttributeError:
        widgets = [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()]
        pb = progressbar.ProgressBar(widgets=widgets, maxval=num_tests).start()
        pb.update(i)
        return pb


if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [pjoin(derec_home, 'fomostos')]
    store_dirs.append(pjoin('/','scratch', 'local1', 'marius', 'doctar_inversion', 'gfdb'))
    store_dirs.append(pjoin('/','scratch', 'local1', 'marius'))

    noisedir = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_noise',
            'displacement')
    time_string = '%s-%s-%s'%time.gmtime()[3:6]
    note = 'with_noise'
    false_store_id = None#'false_doctar_mainland_20Hz'

    file_name = 'robust_check%s_%s.txt'%(time_string, note)
    num_stations = 10
    dz = 2*km
    #num_depths = 11
    num_tests = 1001
    
    engine = LocalEngine(store_superdirs=store_dirs)
    test_type = 'doctar'
    pb = None

    print 'dont add noise'
    add_noise = False
    verbose = True
    debug = False
    write_depth = True

    if test_type=='doctar':
        stf = [[0.,0.2], [0.,1.]]
        #store_id = 'crust2_m5_10Hz'
        # Das ist der neu gebaute store:
        store_id = 'doctar_mainland_20Hz_200m'
        data_base_dir = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_2011-11-01')
        stations_file = 'stations.txt'
        event_file = 'doctar_2011-11-01_quakefile.dat'

    elif test_type=='castor':
        stf = [[0.,1.], [0.,1.]]
        store_id = 'castor'
        data_base_dir = pjoin(derec_home, 'mseeds', 'castor')
        stations_file = 'stations.txt'
        event_file = 'castor_event_2013-10-01.dat'

    phase_ids_start = ['p','P','pP']#'|'.join(du.get_tabulated_phases(engine,
                      #                                 store_id, 
                      #                                 ['p','P']))

    stations = model.load_stations(pjoin(data_base_dir, stations_file))

    targets = du.stations2targets(stations, store_id)

    event = model.Event(load=pjoin(data_base_dir, event_file))
    _ref_source = DCSource.from_pyrocko_event(event)

    if test_type=='doctar':
        targets = filter(lambda x: x.distance_to(_ref_source)<50000., targets)

    #depths = num.linspace(_ref_source.depth-dz, _ref_source.depth+dz, num_depths)
    offset = 3000

    depths = du.drange(_ref_source.depth-offset, _ref_source.depth+offset, 1000)
    #depths=[_ref_source.depth]
    print depths

    ref_source_moment_tensor = _ref_source.pyrocko_moment_tensor()
    location_test_sources_lists = du.make_lots_of_test_events(_ref_source, depths, 
            {'strike':5., 'dip':5., 'rake':5., 'north_shift':1000,
                'east_shift': 1000.}, 
            num_tests,
            func='normal') 
    i=0

    # setup the misfit setup:
    norm = 2
    #taper = trace.CosFader(xfrac=0.333) 
    taper = trace.CosFader(xfade=1.0) 
    
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

    #ok:
    test_case_setup = TestCaseSetup(reference_source=_ref_source,
                                   sources=location_test_sources_lists[0],
                                   targets=targets,
                                   engine=engine, 
                                   store_id=store_id,
                                   misfit_setup=misfit_setup,
                                   source_time_function=stf,
                                   number_of_time_shifts=41,
                                   percentage_of_shift=20.,
                                   phase_ids_start=phase_ids_start,
                                   static_length=1.8,
                                   marker_perc_length=50.0,
                                   marker_shift_frac=0.50,
                                   depths=depths) 
    
    extended_ref_marker, phase_cache = du.chop_ranges(_ref_source, 
                                        targets, 
                                        engine.get_store(store_id),
                                        phase_ids_start,
                                        return_cache=True,
                                        perc=test_case_setup.marker_perc_length,
                                        static_length=test_case_setup.static_length,
                                        t_shift_frac=test_case_setup.marker_shift_frac,
                                        use_cake=True)

    if add_noise:
        noise_fns = glob.glob(noisedir+'/*')
        noise = []
        for fn in noise_fns:
            noise.extend(io.load(fn))
        if not noise:
            print 'wanted to add noise, but didnt find some'
    else:
        noise = None

    reference_seismograms = core.make_reference_trace(_ref_source,
                                                    targets, engine,
                                                    stf,
                                                    noise=noise)

    results = []

    if false_store_id:
        test_case_setup.store_id = false_store_id
        test_case_setup.engine.store_id = false_store_id
        for t in test_case_setup.targets:
            t.store_id = false_store_id

    for location_test_sources in location_test_sources_lists:
        i+=1

        test_case_setup.sources = location_test_sources

        test_case = core.TestCase( test_case_setup )
        test_case.pre_highpass = (2.,0.4)
        test_case.yaml_dump_setup('doctar_setup.yaml')
        test_case.phase_cache = phase_cache
        test_case.set_raw_references(reference_seismograms)

        test_case.set_reference_markers(extended_ref_marker)

        test_case.process()

        best_source, best_misfit = test_case.best_source_misfit()
        angle_diff = best_source.pyrocko_moment_tensor().\
                    angle(ref_source_moment_tensor)

        angle_sign = 1.
        if best_source.strike<=_ref_source.strike:
            angle_sign = -1.

        angle_diff *= angle_sign
        try:
            lateral_shift = num.sqrt(best_source.north_shift**2+best_source.east_shift**2)
            lateral_shift *= best_source.east_shift/abs(best_source.east_shift)
        except ZeroDivisionError:
            print 'both sourcess the same'
            angle_diff = 0
            lateral_shift = 0

        if abs(best_source.depth-_ref_source.depth)<=200:
            if verbose:
                print 'got it, a: %s, north shift: %s'%(angle_diff, lateral_shift)
            if write_depth:
                got_it = best_source.depth
            else:
                got_it = 1
        else:
            if verbose: 
                print 'failed, a: %s, north shift: %s'%(angle_diff, lateral_shift)
            if write_depth:
                got_it = best_source.depth
            else:
                got_it = 0

        if debug:
            op = optics.OpticBase(test_case)
            op.stack_plot()
            plt.show()

        pb = pbar(i, num_tests, pb)

        results.append([lateral_shift, angle_diff, best_misfit, got_it])
        
        if len(results)==10:
            try:
                f = open(file_name, 'a+')
                for line in results:
                    f.write('%s %s %s %s\n'%tuple(line))
            except:
                pass
            finally:
                f.close()
                results = []




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


