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
import socket


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
    if not socket.gethostname()=='Mariuss-MacBook.local':
        store_dirs.append(pjoin('/','scratch', 'local1', 'marius', 
            'doctar_inversion', 'gfdb'))
        store_dirs.append(pjoin('/','scratch', 'local1', 'marius'))

    noisedir = pjoin(derec_home, 'mseeds', 'iris_data', 're-checked_noise')
    time_string = '%s-%s-%s'%time.gmtime()[3:6]
    # MACHT PROBLEME:
    #note = 'false_castor3'
    #note = 'false_castor2'

    note = 'noise_scaled_falsemwp.3'
    file_name = 'robust_check%s_%s.txt'%(time_string, note)
    num_tests = 2000
    use_cake = True
    
    engine = LocalEngine(store_superdirs=store_dirs)
    test_type = 'castor'
    pb = None
    add_noise = True
    verbose = False
    debug = False
    write_depth = True
    false_store_id = None #'false_castor2'
    do_scale = True
    false_magnitude = 0.3

    if test_type=='doctar':
        stf = [[0.,0.1], [0.,1.]]
        store_id = 'doctar_mainland_20Hz'
        data_base_dir = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_2011-11-01')
        stations_file = 'stations.txt'
        event_file = 'doctar_2011-11-01_quakefile.dat'
        phase_ids_start = ['p', 'P'] 

    elif test_type=='castor':
        stf = [[0.,1.], [0.,1.]]
        store_id = 'castor_20Hz'
        data_base_dir = pjoin(derec_home, 'mseeds', 'castor')
        stations_file = 'stations.txt'
        event_file = 'castor_event_2013-10-01.dat'
        phase_ids_start = ['p', 'P', 'Pv2.5p', 'Pv12.5p', 'Pv18.5p', 'Pv20p'] 

    stations = model.load_stations(pjoin(data_base_dir, stations_file))

    targets = du.stations2targets(stations, store_id)

    event = model.Event(load=pjoin(data_base_dir, event_file))
    _ref_source = DCSource.from_pyrocko_event(event)

    if test_type=='doctar':
        targets = filter(lambda x: x.distance_to(_ref_source)<50000., targets)
    
    depths = du.drange(600., 5000., 200)
    print depths
    smaller_magnitude_source = du.clone(_ref_source)

    if false_magnitude:
        smaller_magnitude_source.magnitude = _ref_source.magnitude+false_magnitude
        print 'setting false magnitude to ', smaller_magnitude_source.magnitude

    ref_source_moment_tensor = _ref_source.pyrocko_moment_tensor()
    location_test_sources_lists = du.make_lots_of_test_events(smaller_magnitude_source,
            depths, 
            {('strike','dip', 'rake'):15, ('north_shift', 'east_shift'): 4000}, 
            num_tests,
            func='normal') 
    i=0

    # setup the misfit setup:
    norm = 2
    taper = trace.CosFader(xfrac=0.333) 
    #taper = trace.CosFader(xfade=3.0) 
    
    z, p, k = butter(2, [0.1*num.pi*2, 4.0*num.pi*2.], 
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
                                    number_of_time_shifts=31,
                                    percentage_of_shift=10.,
                                    phase_ids_start=phase_ids_start,
                                    static_length=9.,
                                    marker_perc_length=5.,
                                    marker_shift_frac=0.333,
                                    depths=depths) 

    extended_ref_marker = du.chop_ranges(_ref_source, 
                                        targets, 
                                        engine.get_store(store_id),
                                        phase_ids_start,
                                        perc=test_case_setup.marker_perc_length,
                                        static_length=test_case_setup.static_length,
                                        t_shift_frac=test_case_setup.marker_shift_frac,
                                        use_cake=use_cake)

    if add_noise:
        noise_fns = glob.glob(noisedir+'/*')
        noise = []
        for fn in noise_fns:
            noise.extend(io.load(fn))
        if not noise:
            print 'wanted to add noise, but didnt find some'
    else:
        noise = None

    results = []

    if false_store_id:
        test_case_setup.store_id = false_store_id
        test_case_setup.engine.store_id = false_store_id
        for t in test_case_setup.targets:
            t.store_id = false_store_id

    for location_test_sources in location_test_sources_lists:
        reference_seismograms = core.make_reference_trace(_ref_source,
                                                        targets, engine,
                                                        stf,
                                                        noise=noise)
        i+=1

        test_case_setup.sources = location_test_sources

        test_case = core.TestCase( test_case_setup )
        test_case.set_raw_references(reference_seismograms)

        test_case.set_reference_markers(extended_ref_marker)
        if do_scale:
            test_case.scaling_factors=[1.]
        else:
            test_case.scaling_factors = num.linspace(0.1,2.1,40)

        try:
            test_case.process(verbose=verbose, use_cake=use_cake)
        except meta.OutOfBounds:
            print 'OUT OF BOUNDS... continue'
            continue 

        best_source, best_misfit = test_case.best_source_misfit()
        angle_diff = best_source.pyrocko_moment_tensor().\
                    angle(ref_source_moment_tensor)

        angle_sign = 1.
        if best_source.strike<=_ref_source.strike:
            angle_sign = -1.

        angle_diff *= angle_sign
        lateral_shift = num.sqrt(best_source.north_shift**2+best_source.east_shift**2)
        print best_source.east_shift
        print lateral_shift

        try:
            lateral_shift *= best_source.east_shift/abs(best_source.east_shift)
        except ZeroDivisionError:
            print best_source.east_shift

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
            plt.figure()
            op = optics.OpticBase(test_case)
            op.stack_plot()
            plt.figure()
            op.stack_plot(scaling=test_case.scaling, force_update=True)
            misfit_fig = plt.figure()
            misfit_ax1 = misfit_fig.add_subplot(212)
            misfit_ax1.set_title('scaled')
            op.plot_scaled_misfits(ax=misfit_ax1, 
                                   marker='o', 
                                   color='b', 
                                   lw=0)

            misfit_ax2 = misfit_fig.add_subplot(211)
            misfit_ax2.set_title('un-scaled')
            op.plot_misfits(ax=misfit_ax2, 
                            marker='o', 
                            color='r',
                           lw=0)
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


