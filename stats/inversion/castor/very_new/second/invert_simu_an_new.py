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
import pdb
from derec.inverter import FocusMonteCarlo

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
    num_random_events = 1000

    for nn in range(num_random_events):

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

        note = ''
        file_name = 'very_new_simul%s_%s.txt'%(time_string, note)
        use_cake = False
         
        engine = LocalEngine(store_superdirs=store_dirs)
        test_type = 'castor'
        pb = None
        add_noise = False
        verbose = False
        debug = False
        light_debug = True
        debug_inversion = False
        write_depth = True
        false_store_id = None #'false_castor2'
        do_scale = True
        do_individual_scaling = True
        false_magnitude = None#0.3

        
        #depths = du.drange(1000., 5000., 1000)
        depths = [ 2000.]
        #fine_gridded_depths = du.drange(600., 5000., 200)
        fine_gridded_depths = du.drange(1000., 5000., 1000)

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
        rand_DC = moment_tensor.MomentTensor.random_dc()
 
        _ref_source = DCSource.from_pyrocko_event(event)

        if test_type=='doctar':
            targets = filter(lambda x: x.distance_to(_ref_source)<50000., targets)
        # setup the misfit setup:
        norm = 2
        taper = trace.CosFader(xfrac=0.333) 
        #taper = trace.CosFader(xfade=3.0) 
        
        #z, p, k = butter(2, [0.08*num.pi*2, 4.0*num.pi*2.], 
        z, p, k = butter(2, [0.05*num.pi*2, 2.0*num.pi*2.], 
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
        
        
        du.randomize_DCSource(_ref_source, inplace=True)
        refsourcestr = '0.0 %s %s %s\n'%(_ref_source.strike, _ref_source.dip,\
                _ref_source.rake)
        print _ref_source

        reference_seismograms = core.make_reference_trace(_ref_source,
                                                          targets, 
                                                          engine,
                                                          stf)

        du.randomize_DCSource(_ref_source, inplace=True)
        second_randomized = '%s %s %s\n'%(_ref_source.strike, _ref_source.dip,\
                _ref_source.rake)

        test_case_setup = TestCaseSetup(reference_source=_ref_source,
                                        targets=targets,
                                        engine=engine, 
                                        store_id=store_id,
                                        misfit_setup=misfit_setup,
                                        source_time_function=stf,
                                        number_of_time_shifts=41,
                                        time_shift=0.1,
                                        phase_ids_start=phase_ids_start,
                                        static_length=7.,
                                        marker_perc_length=20.,
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

        tstart = time.time() 
        
        test_case = core.TestCase( test_case_setup )
        test_case.set_raw_references(reference_seismograms)
        test_case.set_reference_markers(extended_ref_marker)
        fmc = FocusMonteCarlo(test_case, debug=debug)
        fmc.steps=3
        fmc.run(1, verbose=True)
        results = fmc.get_results_unformatted()
        elapsed = time.time()-tstart
        f = open(file_name, 'a+')
        f.write('Elapsed time: %s\n' %elapsed)
        f.write(refsourcestr)
        f.write(results)
        f.close()




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


