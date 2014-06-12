from derec.yaml_derec import *
from derec.core import *
from derec.optics import OpticBase
from pyrocko.gf import *
from pyrocko import model, trace, io, gui_util
from scipy.signal import butter

import matplotlib.pyplot as plt
import os
import derec.derec_utils as du
import numpy as num
import glob
pjoin = os.path.join

if __name__ ==  "__main__":

    use_markers =True  
    use_cake = True

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [pjoin(derec_home, 'fomostos')]

    engine = LocalEngine(store_superdirs=store_dirs)

    store_id = 'castor_20Hz'
    stations = model.load_stations(pjoin(derec_home, 'mseeds', 'castor',
                    'stations.txt'))
    event = model.Event(load=pjoin(derec_home, 'mseeds', 'castor',
                    'castor_event_2013-10-01.dat'))
    files = pjoin(derec_home, 'mseeds', 'castor', '2013-10-01T03-32-45',
                                                  'displacement2.mseed')
    if use_markers:
        marker_fn = pjoin(derec_home, 'mseeds','castor','2013-10-01T03-32-45','markers.txt')
        markers = gui_util.Marker.load_markers(marker_fn)

    traces = []
    traces = io.load(files)  #traces.extend(io.load(fn) for fn in glob.glob(files))
    #traces = du.flatten_list(traces)

    map(lambda x: x.highpass(2, 0.05), traces)
    #channels = ['HHE', 'HHN', 'HHZ']

    phase_ids_start = ['p','P']

    targets = du.stations2targets(stations, store_id )
    #targets = du.stations2targets(stations, store_id, channels=channels)
    ref_source = DCSource.from_pyrocko_event(event)

    model = du.get_earthmodel_from_engine(engine, store_id) 

    depths=num.linspace(ref_source.depth-1600, ref_source.depth+1600, 21)

    # Das kann mit als Funktion in TestCaseSetup...
    location_test_sources = du.test_event_generator(ref_source, depths)

    map(lambda x: x.regularize(), location_test_sources)

    rise_time=0.4
    
    stf = [[0.,rise_time],[0.,1.]]

    reference_seismograms = du.make_traces_dict(ref_source, targets, traces)
    targets = reference_seismograms.values()[0].keys()

    # setup the misfit setup:
    norm = 2
    #taper = trace.CosFader(xfrac=0.2) 
    taper = trace.CosFader(xfade=2.2) 
    
    z, p, k = butter(2, [0.5*num.pi*2, 3.0*num.pi*2.], 
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

    test_case_setup = TestCaseSetup(reference_source=ref_source,
                                    sources=location_test_sources,
                                    targets=targets,
                                    engine=engine, 
                                    store_id=store_id,
                                    misfit_setup=misfit_setup,
                                    source_time_function=stf,
                                    number_of_time_shifts=200,
                                    percentage_of_shift=60.,
                                    phase_ids_start=phase_ids_start,
                                    static_length=6.,
                                    marker_perc_length=5.,
                                    marker_shift_frac=0.4,
                                    depths=depths) 

    test_case = TestCase( test_case_setup )
    test_case.scale_minmax = True
    #test_case.blacklist = (('Y7','L004','','HHN'),)

    test_case.set_raw_references(reference_seismograms)

    #markers_dict = du.make_markers_dict(ref_source, targets, markers)

    extended_ref_marker, phase_cache = du.chop_ranges(ref_source, 
                                    targets, 
                                    test_case.store,
                                    phase_ids_start,
                                    picked_phases={},
                                    perc=test_case_setup.marker_perc_length,
                                    static_length=test_case_setup.static_length,
                                    t_shift_frac=test_case_setup.marker_shift_frac,
                                    return_cache=True,
                                    use_cake=use_cake)

    test_case.phase_cache = phase_cache
    test_case.set_reference_markers(extended_ref_marker)

    test_case.process(verbose=True, debug=False)

    ob = OpticBase(test_case)
    ob.stack_plot()
    ob.plot_misfits()
    plt.show()
    
