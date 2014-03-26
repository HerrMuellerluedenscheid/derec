from pyrocko.gf import *
from pyrocko import model, gui_util, pile, trace, moment_tensor, io
from collections import defaultdict
from matplotlib import cm
from matplotlib.mlab import griddata
from gmtpy import griddata_auto
from scipy.signal import butter
from scipy.ndimage import zoom
from guts import *

import matplotlib.transforms as transforms
import time
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.lines as pltlines
import progressbar
import os
from derec import derec_utils as du
from derec.core import *
import numpy as num
import copy
import pdb

pjoin = os.path.join
km = 1000.

if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [derec_home + '/fomostos']

    engine = LocalEngine(store_superdirs=store_dirs)
    store_id = 'castor'
    # load stations from file:
    stations = model.load_stations(pjoin(selfdir,
                            '../reference_stations_castor_selection.txt'))

    markers = gui_util.Marker.load_markers(pjoin(selfdir,
                                                '../reference_marker_castor.txt'))

    phase_ids_start = '|'.join(du.get_tabulated_phases(engine,
                                                       store_id, 
                                                       ['p','P']))
    
    # load stations from file:
    # Event==================================================
    event = filter(lambda x: isinstance(x, gui_util.EventMarker), markers)
    assert len(event) == 1
    event = event[0].get_event()
    event.magnitude = 4.3
    event.moment_tensor = moment_tensor.MomentTensor(
                                    m=num.array([[0.0, 0.0, 1.0],
                                                 [0.0, 0.0, 0.0],
                                                 [0.0, 0.0, 0.0]]))


    # generate stations from olat, olon:
    if not stations:
        print 'Generating station distribution.'
        stations = du.station_distribution((event.lat,event.lon),
                                       [[10000., 4], [130000., 8]], 
                                       rotate={3000.:45, 130000.:0})

    targets = du.stations2targets(stations, store_id)

    model = du.get_earthmodel_from_engine(engine, store_id) 

    #TESTSOURCES===============================================
    
    zoffset= 0.
    ref_source = du.event2source(event, 'DC', strike=37.3, dip=30, rake=-3)

    depths=[1800, 2000, 2200]
    #depths=num.linspace(ref_source.depth-zoffset, ref_source.depth+zoffset, 1)

    print depths, '<- depths'

    location_test_sources = [DCSource(lat=ref_source.lat,
                           lon=ref_source.lon,
                           depth=depth,
                           time=event.time,
                           strike=ref_source.strike,
                           dip=ref_source.dip,
                           rake=ref_source.rake,
                           magnitude=event.magnitude) for depth in depths]

    map(lambda x: x.regularize(), location_test_sources)

    reference_request = make_reference_trace(ref_source,
                                                 targets, 
                                                 engine)

    reference_seismograms = du.response_to_dict(reference_request)

    # setup the misfit setup:
    norm = 2.
    taper = trace.CosFader(xfrac=0.2) 
    
    z, p, k = butter(4, (2.*num.pi*2. ,0.4*num.pi*2.) , 
                                       'bandpass', 
                                       analog=True, 
                                       output='zpk')

    z = num.array(z, dtype=complex)
    p = num.array(p, dtype=complex)
    k = num.complex(k)
    fresponse = trace.PoleZeroResponse(z,p,k)
    fresponse.regularize()

    misfit_setup = trace.MisfitSetup(norm=norm,
                                     taper=taper,
                                     domain='time_domain',
                                     filter=fresponse)
    test_case_dict = {}

    rise_times = num.linspace(0.,10,5)
    for rise_time in rise_times:

        rise_time = 1.
        stf = [[0,rise_time], [0,1]]

        test_case_setup = TestCaseSetup(reference_source=ref_source,
                                        sources=location_test_sources,
                                        targets=targets,
                                        engine=engine, 
                                        store_id=store_id,
                                        misfit_setup=misfit_setup,
                                        source_time_function=stf,
                                        number_of_time_shifts=9,
                                        percentage_of_shift=10.,
                                        phase_ids_start=phase_ids_start) 

        test_case = TestCase( test_case_setup )

        for tr in TestCase.iter_dict(reference_seismograms, only_values=True):
            du.add_random_noise_to_trace(tr, A=0.00001)

        test_case.set_raw_references(reference_seismograms)

        # considering that these traces are 'real' traces. Thus, 
        # stf needs to be applied to raw traces.
        test_case.raw_references = du.apply_stf(test_case.raw_references, 
                                test_case_setup.source_time_function)

        extended_ref_marker = du.chop_ranges(ref_source, 
                                            targets, 
                                            test_case.store,
                                            phase_ids_start,
                                            perc=1.0,
                                            t_shift_frac=0.3)

        test_case.set_reference_markers(extended_ref_marker)

        D = Doer(test_case)
        test_case_dict[rise_time] = test_case


