from collections import defaultdict
import matplotlib.lines as pltlines
import matplotlib.pyplot as plt
import time
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

    i=0
    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [pjoin(derec_home, 'fomostos')]
    if not socket.gethostname()==('love' or 'Mariuss-MacBook.local'):
        store_dirs.append(pjoin('/',
                                'scratch', 
                                'local1', 
                                'marius', 
                                'doctar_inversion', 
                                'gfdb'))

    if not socket.gethostname()=='Mariuss-MacBook.local':
        store_dirs.append(pjoin('/',
                                'scratch', 
                                'local1', 
                                'marius'))

    noisedir = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_noise',
            'displacement')
    time_string = '%s-%s-%s'%time.gmtime()[3:6]
    note = 'no_noise_scaled'
    false_store_id = None#'false_crust1_20Hz_200m'
    false_magnitude = 0.
    do_scale = True 
    individual_scaling = False

    file_name = 'robust_%s_%s.txt'%(time_string, note)
    num_stations = 10
    dz = 2*km
    #num_depths = 11
    num_tests = 2000
    
    engine = LocalEngine(store_superdirs=store_dirs)
    test_type = 'doctar'
    pb = None

    add_noise = False
    verbose = False
    debug = True
    write_depth = True
    write_misfit = True
    write_scaling = True
    check_locations = False
    sdr_range = 0.0001
    lateral_range=0.0001

    if test_type=='doctar':
        stf = [[0.,0.15], [0.,1.]]
        store_id = 'doctar_crust_20Hz_200m'
        #store_id = 'crust2_m5_10Hz'
        # Das ist der neu gebaute store:
        data_base_dir = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_2011-11-01')
        stations_file = 'stations.txt'
        event_file = 'doctar_2011-11-01_quakefile.dat'

    elif test_type=='castor':
        stf = [[0.,1.], [0.,1.]]
        store_id = 'castor'
        data_base_dir = pjoin(derec_home, 'mseeds', 'castor')
        stations_file = 'stations.txt'
        event_file = 'castor_event_2013-10-01.dat'

    phase_ids_start = ['begin']

    stations = model.load_stations(pjoin(data_base_dir, stations_file))

    targets = du.stations2targets(stations, store_id)
    event = model.Event(load=pjoin(data_base_dir, event_file))
    st, di, ra=event.moment_tensor.both_strike_dip_rake()[0]
    _ref_source = DCSource.from_pyrocko_event(event)

    print 'reference source  magnitude: ', _ref_source.magnitude

    if test_type=='doctar':
        targets = filter(lambda x: x.distance_to(_ref_source)<55000., targets)

    ignore = ['L002', 'L003', 'L012','L005', 'L006', 'L007', 'L008']
    targets = filter(lambda x: x.codes[1] not in ignore, targets)
    offset = 4000

    depths = du.drange(1000, 8000, 2000)
    #depths = [5000]
    print 'depths: ', depths
    
    smaller_magnitude_source = du.clone(_ref_source)
    if false_magnitude:
        smaller_magnitude_source.magnitude = _ref_source.magnitude+false_magnitude
        print 'setting false magnitude to ', smaller_magnitude_source.magnitude

    location_test_sources_lists = du.make_lots_of_test_events(smaller_magnitude_source, depths, 
            {('strike', 'dip', 'rake'):sdr_range, ('north_shift', 'east_shift'):
                lateral_range}, 
            num_tests,
            func='uniform') 
    i=0

    if check_locations:
        optics.check_locations(location_test_sources_lists, _ref_source )

    # setup the misfit setup:
    norm = 2
    taper = trace.CosFader(xfrac=0.333) 
    #taper = trace.CosFader(xfade=1.) 
    
    z, p, k = butter(4, [0.5*num.pi*2, 4.0*num.pi*2.], 
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
                                   number_of_time_shifts=0,
                                   time_shift=0.1,
                                   #percentage_of_shift=15.,
                                   phase_ids_start=phase_ids_start,
                                   static_length=3.5,
                                   marker_perc_length=5.0,
                                   marker_shift_frac=0.333,
                                   depths=depths) 
    
    extended_ref_marker, phase_cache = du.chop_ranges(_ref_source, 
                                        targets, 
                                        engine.get_store(store_id),
                                        phase_ids_start,
                                        return_cache=True,
                                        perc=test_case_setup.marker_perc_length,
                                        static_length=test_case_setup.static_length,
                                        t_shift_frac=test_case_setup.marker_shift_frac,
                                        use_cake=False)

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

        
     
    print 'Checked: avoiding immutable seismograms problem...'
    #for location_test_sources in location_test_sources_lists:
    for i in range(num_tests):
        du.randomize_DCSource(_ref_source, inplace=True)
        ref_source_moment_tensor = _ref_source.pyrocko_moment_tensor()
        if false_store_id:
            for t in test_case_setup.targets:
                t.store_id = store_id 
        reference_seismograms = core.make_reference_trace(_ref_source,
                                                        targets, engine,
                                                        stf,
                                                        noise_scale=1.,
                                                        noise_type='natural',
                                                        noise=noise)

        if false_store_id:
            test_case_setup.store_id = false_store_id
            test_case_setup.engine.store_id = false_store_id
            for t in test_case_setup.targets:
                t.store_id = false_store_id
        _ref_source.regularize()
        print _ref_source
        du.invert_DC(_ref_source, inplace=True)
        _ref_source.regularize()
        print _ref_source
        location_test_sources = du.make_lots_of_test_events(_ref_source, depths, 
                {('strike', 'dip', 'rake'):sdr_range, ('north_shift',
                    'east_shift'): lateral_range}, 
                1,
                func='uniform')[0]
        i+=1

        test_case_setup.sources = location_test_sources

        test_case = core.TestCase( test_case_setup )
        test_case.pre_highpass = (2.,0.5)
        test_case.phase_cache = phase_cache
        
        if do_scale:
            test_case.scaling_factors = num.linspace(0.1,2.1, 40)
        else:
            test_case.scaling_factors = [1.]

        test_case.individual_scaling = individual_scaling
        test_case.reduce_half_rise = False

        test_case.set_raw_references(reference_seismograms)

        test_case.set_reference_markers(extended_ref_marker)

        try:
            test_case.process(verbose=verbose)
        except meta.OutOfBounds:
            continue

        best_source, best_misfit = test_case.best_source_misfit(verbose=verbose)
        angle_diff = best_source.pyrocko_moment_tensor().\
                    angle(ref_source_moment_tensor)

        
        best_scaling = test_case.scaling[best_source]
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
            print 'best depth: ', best_source.depth
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
        if write_misfit:
            data2write = [lateral_shift, angle_diff, best_misfit, got_it,
                    best_scaling]
        else:
            data2write = [lateral_shift, angle_diff, best_misfit, got_it]

        results.append(data2write)
        
        if len(results)==5:
            try:
                f = open(file_name, 'a+')
                for line in results:
                    if write_misfit and write_scaling:
                        f.write('%s %s %s %s %s\n'%tuple(line))
                    else:
                        f.write('%s %s %s %s\n'%tuple(line))
            except:
                raise
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


