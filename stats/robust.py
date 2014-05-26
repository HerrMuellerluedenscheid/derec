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
    #store_dirs.append(pjoin('/','scratch', 'local1', 'marius', 'doctar_inversion', 'gfdb'))
    noisedir = pjoin(derec_home, 'mseeds', 'iris_data', 're-checked_noise')
    time_string = '%s-%s-%s'%time.gmtime()[3:6]
    file_name = 'robust_check_results_%s.txt'%time_string
    num_stations = 10
    stf = [[0.,1.], [0.,1.]]
    dz = 2*km
    num_depths = 9
    num_tests = 10
    
    engine = LocalEngine(store_superdirs=store_dirs)
    test_type = 'doctar'
    pb = None
    add_noise = True
    verbose = True
    
    if test_type=='doctar':
        store_id = 'doctar_mainland_20Hz'
        data_base_dir = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_2011-11-01')
        stations_file = 'stations.txt'
        event_file = 'doctar_2011-11-01_quakefile.dat'
        phase_ids_start = '|'.join(du.get_tabulated_phases(engine,
                                                           store_id, 
                                                           ['p','P']))

    elif test_type=='castor':
        store_id = 'castor'
        data_base_dir = pjoin(derec_home, 'mseeds', 'castor')
        stations_file = 'reference_stations_castor.txt'
        event_file = 'castor_event_2013-10-01.dat'
        phase_ids_start = '|'.join(du.get_tabulated_phases(engine,
                                                           store_id, 
                                                           ['p','P']))

    stations = model.load_stations(pjoin(data_base_dir, stations_file))

    targets = du.stations2targets(stations, store_id)
    event = model.Event(load=pjoin(data_base_dir, event_file))
    _ref_source = DCSource.from_pyrocko_event(event)
    depths = num.linspace(_ref_source.depth-dz, _ref_source.depth+dz, num_depths)
    ref_source_moment_tensor = _ref_source.pyrocko_moment_tensor()
    location_test_sources_lists = du.make_lots_of_test_events(_ref_source, depths, 
            {'strike':10, 'dip':10, 'rake':10, 'north_shift':5000 }, 
            num_tests,
            func='normal') 
    i=0

    # setup the misfit setup:
    norm = 2
    taper = trace.CosFader(xfrac=0.2) 
    
    z, p, k = butter(2, [0.001*num.pi*2, 2.0*num.pi*2.], 
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

    test_case_setup = TestCaseSetup(reference_source=_ref_source,
                                    sources=location_test_sources_lists[0],
                                    targets=targets,
                                    engine=engine, 
                                    store_id=store_id,
                                    misfit_setup=misfit_setup,
                                    source_time_function=stf,
                                    number_of_time_shifts=9,
                                    percentage_of_shift=10.,
                                    phase_ids_start=phase_ids_start,
                                    static_length=4.,
                                    marker_perc_length=50.,
                                    marker_shift_frac=0.2,
                                    depths=depths) 

    extended_ref_marker = du.chop_ranges(_ref_source, 
                                        targets, 
                                        engine.get_store(store_id),
                                        phase_ids_start,
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

    reference_seismograms.values()[0].values()[0].snuffle()


    results = []

    for location_test_sources in location_test_sources_lists:
        i+=1

        test_case_setup.sources = location_test_sources

        test_case = core.TestCase( test_case_setup )
        test_case.set_raw_references(reference_seismograms)

        test_case.set_reference_markers(extended_ref_marker)

        test_case.process()
        

        best_source, best_misfit = test_case.best_source_misfit()
        angle_diff = best_source.pyrocko_moment_tensor().\
                    angle(ref_source_moment_tensor)

        angle_sign = 1.
        if best_source.pyrocko_moment_tensor().strike<=ref_source_moment_tensor.strike:
            angle_sign = -1.

        angle_diff *= angle_sign
        lateral_shift = num.sqrt(best_source.north_shift**2+best_source.east_shift**2)
        lateral_shift *= best_source.east_shift/abs(best_source.east_shift)

        if abs(best_source.depth-_ref_source.depth)<=200:
            if verbose:
                print 'got it, a: %s, north shift: %s'%(angle_diff, lateral_shift)
            got_it = 1
        else:
            if verbose: 
                print 'failed, a: %s, north shift: %s'%(angle_diff, lateral_shift)
            got_it = 0

        if verbose:
            op = optics.OpticBase(test_case)
            op.stack_plot()
            plt.show()

        if not verbose:
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


