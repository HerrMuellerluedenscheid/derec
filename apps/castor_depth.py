from derec.yaml_derec import *
from derec.core import *
from derec.optics import OpticBase, plot_misfit_dict
from derec.inverter import FocusMonteCarlo
from pyrocko.gf import *
from pyrocko import model, trace, io, gui_util, guts
from scipy.signal import butter

import matplotlib.pyplot as plt
import os
import derec.derec_utils as du
import numpy as num
import glob
import socket
pjoin = os.path.join

if __name__ ==  "__main__":

    use_markers = True  
    use_cake = True
    test_type  = 'castor'
    invert =False

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

    if not socket.gethostname() in ('Mariuss-MacBook.local', 'M.local'):
        store_dirs.append(pjoin('/',
                                'scratch',
                                'local1',
                                'marius'))
        store_dirs.append('/data/share/u253/doctar/gfdb_stores')
    engine = LocalEngine(store_superdirs=store_dirs)
 

    if test_type == 'castor':
        store_id = 'castor_20Hz'
        stations = model.load_stations(pjoin(derec_home, 'mseeds', 'castor',
                        'stations_selection.txt'))
        event = model.Event(load=pjoin(derec_home, 'mseeds', 'castor',
                        'castor_event_2013-10-01.dat'))
        files = [pjoin(derec_home, 'mseeds', 'castor', '2013-10-01T03-32-45',
                                                      'displacement2.mseed')]
        #test_case_setup = guts.load(filename=derec_home+'/stats/castor_results/castor_standard.yaml') 
        test_case_setup = guts.load(filename=derec_home+'/stats/castor_results/castor_mod.yaml') 
        depths=num.arange(1000., 5000.,200.)

        if use_markers:
            marker_fn = pjoin(derec_home, 
                              'mseeds',
                              'castor',
                              '2013-10-01T03-32-45',
                              'markers.txt')
            markers = gui_util.Marker.load_markers(marker_fn)

    elif test_type == 'doctar':
        print 'ACHTUNG! DAS IST NICHT DAS CRUST MODEL SONDERN DAS MAINLAND'
        store_id = 'doctar_crust_20Hz_200m'
        data_base_dir = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_2011-11-01')
        stations = model.load_stations(data_base_dir+'/stations.txt')
        event_file = 'doctar_2011-11-01_quakefile.dat'
        files = glob.glob(pjoin(data_base_dir, 'restituted')+'/*')
        event = model.Event(load=pjoin(data_base_dir, 
                                       'doctar_2011-11-01_quakefile.dat'))
        test_case_setup = guts.load(filename=pjoin(derec_home, 
            'stats/doctar_results/crust/random_mts'+\
                                                    '/doctar_standard.yaml')) 
        depths=num.arange(1000., 8200.,200.)
        if use_markers:
            print 'MARKERS EINTRAGE!!'
            marker_fn = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_2011-11-01',
                              'doctar_markers_111101.txt')
            markers = gui_util.Marker.load_markers(marker_fn)

    traces = [io.load(f) for f in files]
    traces = du.flatten_list(traces)

    #map(lambda x: x.highpass(2, 0.05), traces)
    if test_type=='doctar':
        print 'highpassing traces'
        map(lambda x: x.highpass(2, 0.5), traces)
    else:
        print 'highpassing traces'
        map(lambda x: x.highpass(4, 0.1), traces)

    print 'event magnitude: ', event.magnitude

    targets = du.stations2targets(stations, store_id )
    ref_source = DCSource.from_pyrocko_event(event)
    print 'ref source magnitude: ', ref_source.magnitude
    if test_type=='castor':
        ref_source.magnitude*=0.65
    model = du.get_earthmodel_from_engine(engine, store_id) 
    print 'using depths %s'%depths

    # Das kann mit als Funktion in TestCaseSetup...
    location_test_sources = du.test_event_generator(ref_source, depths)
    map(lambda x: x.regularize(), location_test_sources)

    norm = 2
    taper = trace.CosFader(xfrac=0.333)
    #taper = trace.CosFader(xfade=3.0) 

    z, p, k = butter(2, [0.05*num.pi*2, 3.0*num.pi*2.],
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
    
    
    print ' NACHSCHAUEN OB DAS WIRKLICH DAS RICHTIGE CASTOR STANDARD FILE IST!'
    test_case_setup.targets = targets
    test_case_setup.engine = engine
    test_case_setup.reference_source = ref_source
    test_case_setup.store_id = store_id
    #test_case_setup.number_of_time_shifts = 0
    test_case_setup.percentage_of_shift= 5.
    test_case_setup.validate()
    test_case_setup.misfit_setup = misfit_setup
    print 'MISFIT SETUP AUSGETAUSCHT !!!!'
    test_case_setup.sources = location_test_sources
    reference_seismograms = du.make_traces_dict(ref_source, test_case_setup.targets, traces)
    targets = reference_seismograms.values()[0].keys()
    
    test_case = TestCase( test_case_setup )
    test_case.individual_scaling = True
    test_case.scaling_factors = num.arange(0.1, 4., 0.1)
    test_case.blacklist = (('Y7','L004','','HHN'),('ES','EMOS','','N'))

    test_case.set_raw_references(reference_seismograms)

    markers_dict = du.make_markers_dict(ref_source, targets, markers)

    extended_ref_marker, phase_cache = du.chop_ranges(ref_source, 
                                    targets, 
                                    test_case.store,
                                    test_case_setup.phase_ids_start,
                                    picked_phases=markers_dict,
                                    perc=test_case_setup.marker_perc_length,
                                    static_length=test_case_setup.static_length,
                                    t_shift_frac=test_case_setup.marker_shift_frac,
                                    return_cache=True,
                                    use_cake=use_cake)

    test_case.phase_cache = phase_cache
    test_case.set_reference_markers(extended_ref_marker)

    if not invert:
        test_case.process(verbose=True, debug=False, use_cake=use_cake)
    if invert:
        print 'REDUCING NUMBER OF DEPTHS'
        test_case.test_case_setup.depths = [2000.,5000.,7000.,]
        inv = FocusMonteCarlo(test_case)
        #inv.set_focus([360., 90., 20., 5])
        inv.set_focus([360., 90., 20., 5])
        inv.run(1, verbose=True)
        results = inv.get_results_unformatted()
        f = open('inversion.txt', 'w')
        f.write(results)
        f.close()



    ob = OpticBase(test_case)
    plt.figure()
    plot_misfit_dict(test_case.misfits, 
                     test_case.scaled_misfits, 
                     scaling=test_case.scaling)
    plt.figure()
    ob.stack_plot()
    plt.figure()
    ob.stack_plot(scaling=test_case.scaling, force_update=True)
    #ob.plot_misfits()
    plt.show()
    
