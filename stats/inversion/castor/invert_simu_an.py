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



def accept(bestmf, lastbestmf, i, num_inversions, T=None):
    if bestmf<lastbestmf:
        return 1 

    if not T:
        l = 1.5
        T = num.linspace(0,1, num_inversions)
        T = T*num.exp(-l*T)

        e = 3.5
        T_increase = num.linspace(0,1, num_inversions)
        T_increase = T_increase*num.exp(-l*T_increase)

    T_c = num.random.random(1)
    if T_c > T[i]:
        T_i= num.random.random(1)
        if T_i<T_increase[i]:
            # this means increase temperature on step
            return 2
        else:
            return 0 
    else:
        return 1 
    

if __name__ ==  "__main__":

    num_random_events = 200

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
        #note = 'false_castor3'
        #note = 'false_castor2'

        note = 'noise_scaled'
        file_name = 'robust_check%s_%s.txt'%(time_string, note)
        num_tests = 21
        use_cake = False
        
        engine = LocalEngine(store_superdirs=store_dirs)
        test_type = 'castor'
        pb = None
        add_noise = True
        verbose = False
        debug = False
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

        num_inversions = 5
        max_sdr_freedom = 180
        skip_threshold = 0.9
        skip_threshold_num = 5
        circles_before_break = 30

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
        _ref_source.strike, _ref_source.dip, _ref_source.rake = rand_DC.both_strike_dip_rake()[0]

        if test_type=='doctar':
            targets = filter(lambda x: x.distance_to(_ref_source)<50000., targets)
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
        

        last_best = du.clone(_ref_source)
        last_best.strike= 0.
        last_best.dip=90.
        last_best.rake=0.
        last_best.depth = 2000
        init_angle = _ref_source.pyrocko_moment_tensor().angle(last_best.pyrocko_moment_tensor())
        try:
            f = open(file_name, 'a+')
            f.write('%s %s %s %s %s %s\n'%(999, _ref_source.strike,
                                                      _ref_source.dip,
                                                      _ref_source.rake,
                                                      0.,
                                                      _ref_source.depth))
            f.write('%s %s %s %s %s %s\n'%(999, last_best.strike,
                                                      last_best.dip,
                                                      last_best.rake,
                                                      init_angle,
                                                      last_best.depth))
        except:
            raise
        finally:
            f.close()
        print 'ref source is  s %s, d %s, r %s,' %(_ref_source.strike,
                                    _ref_source.dip, _ref_source.rake) 
        print 'input last best source is  s %s, d %s, r %s,' %(last_best.strike,
                                    last_best.dip, last_best.rake) 

        if false_magnitude:
            last_best.magnitude = last_best.magnitude+false_magnitude
            print 'setting false magnitude to ', last_best.magnitude
        
        sdr_freeddom = num.linspace(1, max_sdr_freedom, num_inversions)[::-1]
        i=0
        skipped = 0
        ref_source_moment_tensor = last_best.pyrocko_moment_tensor()
        lastbestmf = 0.9

        best_mf_of_step = defaultdict()
        i_total = 0
        while i<=num_inversions-1:

            if i==num_inversions-1:
                depths = fine_gridded_depths
            print 'step %s of %s ' %(i, num_inversions)
            ref_source = du.clone(last_best)

            print 'new freedom: ', sdr_freeddom[i]
            if i==0:
                isfirst=True
                num_tests_used = num_tests*2
            else:
                num_tests_used = num_tests
                isfirst=False

            location_test_sources_lists = du.make_lots_of_test_events(ref_source,
                    depths, 
                    {('strike','dip', 'rake'):sdr_freeddom[i]}, 
                    num_tests_used,
                    func='uniform', 
                    isfirst=isfirst) 

            if isfirst:
                print 'using the following sdr options:'
                for loct in location_test_sources_lists:
                    print loct[0].strike, loct[0].dip, loct[0].rake
            if debug_inversion:
                optics.check_locations(ref_source, location_test_sources_lists)

            test_case_setup = TestCaseSetup(reference_source=ref_source,
                                            sources=location_test_sources_lists[0],
                                            targets=targets,
                                            engine=engine, 
                                            store_id=store_id,
                                            misfit_setup=misfit_setup,
                                            source_time_function=stf,
                                            number_of_time_shifts=31,
                                            time_shift=0.3,
                                            phase_ids_start=phase_ids_start,
                                            static_length=9.,
                                            marker_perc_length=5.,
                                            marker_shift_frac=0.333,
                                            depths=depths) 

            extended_ref_marker = du.chop_ranges(ref_source, 
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
            best_sources = defaultdict()
            all_the_best = []
            for location_test_sources in location_test_sources_lists:
                reference_seismograms = core.make_reference_trace(ref_source,
                                                                targets, engine,
                                                                stf,
                                                                noise=noise)

                test_case_setup.sources = location_test_sources

                test_case = core.TestCase( test_case_setup )
                test_case.set_raw_references(reference_seismograms)

                test_case.set_reference_markers(extended_ref_marker)
                if do_scale:
                    test_case.scaling_factors=[1.]
                else:
                    test_case.scaling_factors = num.linspace(0.1,2.1,40)

                if do_individual_scaling:
                    test_case.individual_scaling = True

                try:
                    test_case.process(verbose=verbose, use_cake=use_cake)
                except meta.OutOfBounds:
                    print 'OUT OF BOUNDS... continue'
                    continue 

                best_source, best_misfit = test_case.best_source_misfit()


                best_sources[best_misfit] = best_source
                angle_diff = best_source.pyrocko_moment_tensor().\
                            angle(ref_source_moment_tensor)

                #lateral_shift = num.sqrt(best_source.north_shift**2+\
                #        best_source.east_shift**2)
              
                #try:
                #    lateral_shift *= best_source.east_shift/abs(best_source.east_shift)
                #except ZeroDivisionError:
                #    print best_source.east_shift

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

            bestmf = min(best_sources.keys())

            best_angle = best_source.pyrocko_moment_tensor().angle(ref_source_moment_tensor)
            if debug_inversion:
                print 'best source SDR: ', best_source.strike,\
                                         best_source.dip,\
                                        best_source.rake

            do_accept = accept(bestmf, lastbestmf, i, num_inversions)
            if do_accept == 1:
                lastbestmf = bestmf
                best_mf_of_step[i]=bestmf
                i+=1

                if debug_inversion: print 'accepted'
            if do_accept==0: 

                if debug_inversion: print 'not accepted'
                continue
            if do_accept==2: 
                if debug_inversion: print 'not accepted, increase T'
                if i>0:
                    i-=1
                else:
                    i=0
                bestmf = best_mf_of_step[i]
                continue


            last_best = best_sources[bestmf]
            all_the_best.append([bestmf, last_best])
            print all_the_best
                
            try:
                f = open(file_name, 'a+')
                f.write('%s %s %s %s %s %s\n'%(bestmf, last_best.strike,
                                                          last_best.dip,
                                                          last_best.rake,
                                                          best_angle,
                                                          last_best.depth))
            except:
                raise
            finally:
                f.close()

            if i_total==circles_before_break:
                break




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


