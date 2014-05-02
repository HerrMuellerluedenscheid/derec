import numpy as num
import os
import matplotlib.pyplot as plt
import copy
from collections import defaultdict

from pyrocko.gf.meta import ConfigTypeA
from pyrocko.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko.gf import *
from pyrocko import *
from pyrocko.guts import load_string, String, Float, Int, Dict
from pyrocko.guts_array import Array

from derec import derec_utils as du
from derec import core
from derec import optics
from derec.yaml_derec import *

pjoin = os.path.join
fail_message = 'Need to load a setup, first'
derec_home = os.environ["DEREC_HOME"]
store_dirs = [derec_home + '/fomostos']
#engine = LocalEngine(store_superdirs=store_dirs)
#store_id = 'castor'
k = 1000.
misfit = None


class Derec(Snuffling):
    '''
    
    '''
    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Depth Relocation')
        self.add_parameter(Choice('Store: ', 'store_id_choice',
                                    'Need to select a superdir first.',
                                    ['Need to select a superdir first.']))
        self.add_parameter(Param('Static marker length [s]', 'static_length', \
                3, 1, 10))
        self.add_parameter(Param('Relative Marker length [%]',\
                'marker_perc_length', 10., 0., 100.))
        self.add_parameter(Param('Minimal Depth [m]', 'z_min', 1*k, 0., 100*k))
        self.add_parameter(Param('Maximal Depth [m]', 'z_max', 10*k, 0., 100*k))
        self.add_parameter(Param('Number of Depths', 'num_depths', 3, 1, 21))
        self.add_parameter(Param('rel. Shift of Marker []', 'marker_shift_frac',\
                0.3, 0., 1.))
        self.add_parameter(Param('Number of Time Shifts', 'num_time_shifts',\
                6, 0, 30))
        self.add_parameter(Param('rel. trace shift [%]', 'perc_of_shift', 20.,\
                0., 100.))
        self.add_parameter(Param('Rise Time [s]', 'rise_time', 1., 0.5, 5.))
        self.add_parameter(Param('num inversion steps', 'num_inversion_steps',\
                2., 1, 10.))
        self.add_parameter(Switch('Post invert', 'post_invert', True))
        self.add_parameter(Switch('Pre invert', 'pre_invert', True))
        self.add_trigger('Load Misfit Setup', self.load_misfit_setup) 
        self.add_trigger('Load Default Setup', self.load_setup) 
        self.add_trigger('Generate Markers', self.generate_markers) 
        self.add_trigger('Select Store Directory', self.setup_id_choice)
        self.add_trigger('Save result', self.save_result)
        self.add_trigger('Save setup', self.save_setup)
        self.set_live_update(False)

        self.test_case_setup = None
        self.phase_ids_start = 'p|P'
        self._store_ids = []
        self.sources = []
        self.targets = []
        self.reference_source = None

    def call(self):
        '''Main work routine of the snuffling.'''

        self.cleanup()
        if not self.test_case_setup:
            self.fail(fail_message)

        self.active_event, self.stations = self.get_active_event_and_stations()
        if not self.targets:
            self.targets = du.stations2targets(self.stations, \
                    self.store_id_choice)

        if not self.reference_source:
            self.reference_source = du.event2source(self.active_event, 'DC')

        _depths = num.linspace(self.z_min, self.z_max, self.num_depths)

        if not self.pre_invert:
            depths = _depths
        else:
            depths = [self.reference_source.depth]
        sources = du.test_event_generator(self.reference_source, depths)
        traces = self.get_pile().all()
        self.traces_dict = du.make_traces_dict(self.reference_source, \
                self.targets,
                traces)

        stf = [[0., self.rise_time],[0.,1.]]
        # TODO: Qt4 wie directory waehlen fuer engine dirs
        test_case_setup = TestCaseSetup(reference_source=self.reference_source,
                                        sources=sources,
                                        targets=self.targets,
                                        engine=self.engine,
                                        store_id=self.store_id_choice,
                                        misfit_setup=self.misfit_setup,
                                        source_time_function=stf,
                                        number_of_time_shifts=int(\
                                                self.num_time_shifts),
                                        percentage_of_shift=self.perc_of_shift,
                                        phase_ids_start=self.phase_ids_start,
                                        static_length=self.static_length,
                                        marker_perc_length=\
                                                self.marker_perc_length,
                                        marker_shift_frac=\
                                                self.marker_shift_frac,
                                        depths=depths)
        test_case_setup.regularize()
        if self.pre_invert:
            test_case = self.invert(setup=test_case_setup, 
                                    source=self.reference_source,
                                    traces_dict=self.traces_dict,
                                    marker_dict=self.ref_markers_dict)

            test_case.test_case_setup.depths = _depths

            best_source, minmf = test_case.best_source_misfit()
            test_case_setup.reference_source = best_source
            test_case_setup.sources = du.test_event_generator(\
                    best_source, _depths)


        else:
            test_case = core.TestCase(test_case_setup)
            test_case.set_raw_references(self.traces_dict)
            test_case.set_reference_markers(self.ref_markers_dict)
        test_case.process()

        tmp_out_dir = self.tempdir()

        fig = self.figure()
        optics.plot_misfit_dict(test_case.misfits, ax=fig.gca())
        fig.canvas.draw()

        optic = optics.OpticBase(test_case)
        #fig = self.figure()
        optic.stack_plot()
        plt.show()
        #for ax in axs.values():
        #    fig.add_axes(ax)
        #
        #fig.canvas.draw()
    
        self.dumped_results = test_case.yaml_dump(pjoin(tmp_out_dir, \
                'derec_results.yaml'))
        self.dumped_setup = test_case.yaml_dump_setup(pjoin(tmp_out_dir, \
                'derec_setup.yaml'))

        last_misfit = min(test_case.misfits.values())
        misfit = last_misfit
        for k,v in test_case.misfits.iteritems():
            if v==last_misfit:
                depth = k.depth
                break

        if self.post_invert:
            best_source, last_misfit = test_case.best_source_misfit()
            test_case = self.invert(setup=test_case_setup, 
                                    source=best_source,
                                    last_test_case=test_case,
                                    marker_dict=self.ref_markers_dict)

    def invert(self, setup, source, last_test_case=None,\
            traces_dict=None, marker_dict=None):
        """
        *setup* will not be modified but copied.
        *source* is the source to start with
        """
        sdr = dict(zip(['strike', 'dip', 'rake'],
            [source.strike, source.dip, source.rake]))

        dist=1.
        setup = core.TestCaseSetup(**setup.__dict__)
        setup.depths = [source.depth]
        setup.sources = [source]
        if last_test_case:
            s_, last_misfit = last_test_case.best_source_misfit()
        else:
            last_misfit = 0.

        inversion_misfits = defaultdict()

        for i in xrange(self.num_inversion_steps):
            test_case = core.TestCase(setup)
            if traces_dict:
                test_case.set_raw_references(traces_dict)
            if marker_dict:
                test_case.set_reference_markers(marker_dict)
            test_case.process()
            best_source, misfit = test_case.best_source_misfit()
            inversion_misfits[copy.copy(best_source)] = misfit
            if i==self.num_inversion_steps-1 or abs(misfit-last_misfit)<0.001:
                print misfit, last_misfit
                print 'finised inversion after %s steps'%(i+1)
                source.regularize()
                source.validate()
                best_source.regularize()
                best_source.validate()
                print 'ref source: ', source
                print 'best match: ', best_source
                break

            if misfit>last_misfit:
                grads = self.sdr_grad(misfit, setup)
                last_misfit = misfit 
            sdr = self.new_sdr(sdr, grads, dist)

            for k,v in sdr.items():
                 setattr(setup.sources[0], k, v)
            assert len(setup.sources)==1

            test_case.drop_data('raw_candidates')
            test_case.set_setup(setup)
            print test_case.misfits,'<<< all\n'
            assert len(test_case.misfits)==1

        print '----INVERSION RESULTS:'
        for k,v in inversion_misfits.items():
            k.regularize()
            print '\n%s, %s \n'%(k, v)

        return test_case
                

    def sdr_grad(self, origin_misfit, origin_setup):
        save_assert = origin_setup.reference_source.__dict__

        d_deg = 1.
        mfs = []
        for i, attr in enumerate(['strike', 'dip', 'rake']):
            setup = core.TestCaseSetup(**origin_setup.__dict__)
            assert(setup.reference_source.strike==save_assert['strike'])
            setup.reference_source.regularize()
            setup.depths = [setup.reference_source.depth]

            setattr(setup.reference_source, attr,\
                            getattr(setup.reference_source, attr)+d_deg)

            sources = du.test_event_generator(setup.reference_source,\
                            setup.depths)
            setup.sources = sources
            setup.validate()
            test_case = core.TestCase(setup)
            test_case.set_raw_references(self.traces_dict)
            test_case.set_reference_markers(self.ref_markers_dict)
            test_case.process()
            mfs.append(origin_misfit-test_case.misfits.values()[0])

        direction = num.array(mfs)
        direction = direction/num.linalg.norm(direction)#*-1
        print 'new sdr direction: ', direction
        return direction

    def new_sdr(self, origin, grad, distance):
        return dict(zip(origin.keys(), origin.values()+grad*distance))

    def setup_id_choice(self):
        store_ids = self.input_filename(caption='Select a store') 

    def load_misfit_setup(self):
        fn = self.input_filename(caption='Select a misfit setup')
        f = open(fn,'r')
        self.misfit_setup = load_string(f.read())
        f.close()

    def load_setup(self):
        fn = self.input_filename(caption='Select a setup')
        f = open(fn,'r')
        self.test_case_setup = load_string(f.read())
        f.close()
        
        self.set_parameter('static_length', \
                self.test_case_setup.static_length)
        self.set_parameter('num_time_shifts', \
                float(self.test_case_setup.number_of_time_shifts))
        self.set_parameter('num_depths', \
                len(self.test_case_setup.depths))
        self.set_parameter('z_min', min(self.test_case_setup.depths))
        self.set_parameter('z_max', max(self.test_case_setup.depths))
        self.set_parameter('rise_time', \
                self.test_case_setup.source_time_function[0][1])
        self.set_parameter('marker_perc_length', \
                self.test_case_setup.marker_perc_length)
        self.set_parameter('marker_shift_frac', \
                self.test_case_setup.marker_shift_frac)
        self.set_parameter('perc_of_shift',\
                self.test_case_setup.percentage_of_shift)
        self.phase_ids_start = self.test_case_setup.phase_ids_start
        self.engine = self.test_case_setup.engine
        self.misfit_setup = self.test_case_setup.misfit_setup
        self._store_ids.extend(self.test_case_setup.engine.get_store_ids())
        self.set_parameter_choices('store_id_choice', self._store_ids) 

    def save_result(self):
        fn = self.output_filename('Save Results')
        f = open(fn, 'r')
        f.write(self.dumped_results)
        f.close()

    def save_setup(self):
        fn = self.output_filename('Save Setup')
        f = open(fn, 'r')
        f.write(self.dumped_setup)
        f.close()

    def generate_markers(self):
        try:
            self.viewer.remove_markers(self.markers)
        except AttributeError:
            self.viewer = self.get_viewer()

        self.active_event, self.stations = self.get_active_event_and_stations()

        if not self.targets:
            self.targets = du.stations2targets(self.stations, \
                    self.store_id_choice)

        if not self.reference_source:
            self.reference_source = du.event2source(self.active_event, 'DC')
        self.ref_markers_dict = du.chop_ranges(self.reference_source,
                       self.targets,
                       self.engine.get_store(self.store_id_choice),
                       self.phase_ids_start,
                       perc=self.marker_perc_length,
                       static_length=self.static_length,
                       t_shift_frac=self.marker_shift_frac,
                       use_cake=True)

        self.markers = self.ref_markers_dict.values()[0].values()
        self.add_markers(self.markers)

                
def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ Derec() ]


