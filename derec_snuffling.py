import numpy as num
import tempfile
import os
import pyrocko
import matplotlib.pyplot as plt

from pyrocko.gf.meta import ConfigTypeA
from pyrocko.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko import trace
from pyrocko.gf import *
from pyrocko.gf.meta import ConfigTypeA
from pyrocko import *
from pyrocko.guts import load_string, String, Float, Int, Dict
from pyrocko.guts_array import Array

from derec import derec_utils as du
from derec import core
from derec import optics
from derec.yaml_derec import *
from derec import yaml_derec

pjoin = os.path.join
fail_message = 'Need to load a setup, first'
derec_home = os.environ["DEREC_HOME"]
store_dirs = [derec_home + '/fomostos']
engine = LocalEngine(store_superdirs=store_dirs)
store_id = 'castor'
k = 1000.


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

        depths = num.linspace(self.z_min, self.z_max, self.num_depths)
        sources = du.test_event_generator(self.reference_source, depths)
        stf = [[0., self.rise_time],[0.,1.]]
        # TODO: Qt4 wie directory waehlen fuer engine dirs
        test_case_setup = TestCaseSetup(reference_source=self.reference_source,
                                        sources=sources,
                                        targets=self.targets,
                                        engine=self.engine,
                                        store_id=self.store_id_choice,
                                        misfit_setup=self.misfit_setup,
                                        source_time_function=stf,
                                        number_of_time_shifts=int(self.num_time_shifts),
                                        percentage_of_shift=self.perc_of_shift,
                                        phase_ids_start=self.phase_ids_start,
                                        static_length=self.static_length,
                                        marker_perc_length=self.marker_perc_length,
                                        marker_shift_frac=self.marker_shift_frac,
                                        depths=depths)

        test_case = core.TestCase(test_case_setup)
        traces = self.get_pile().all()
        traces_dict = du.make_traces_dict(self.reference_source, self.targets,\
                traces)
        test_case.set_raw_references(traces_dict)
        test_case.set_reference_markers(self.ref_markers_dict)
        core.Doer(test_case)

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
                       self.engine.get_store(store_id),
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


