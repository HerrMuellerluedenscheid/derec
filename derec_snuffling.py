import numpy as num
import tempfile
import os
import pyrocko
from pyrocko.gf.meta import ConfigTypeA
from pyrocko.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko import trace
from pyrocko.gf import *
from pyrocko.gf.meta import ConfigTypeA
from pyrocko import *
from derec import derec_utils as du
from derec import core
from derec.yaml_derec import *
from derec import yaml_derec
from pyrocko.guts import load_string, String, Float, Int, Dict
from pyrocko.guts_array import Array

pjoin = os.path.join
fail_message = 'Need to load a setup, first'
guts_prefix = 'pyrocko.gf.meta'
derec_home = os.environ["DEREC_HOME"]
store_dirs = [derec_home + '/fomostos']
engine = LocalEngine(store_superdirs=store_dirs)
store_id = 'castor'
engine.regularize()
print engine.get_store(store_id)


class Derec(Snuffling):
    '''
    '''
    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Depth Relocation')
        self.add_parameter(Choice('Store: ', 'store_id',
                                    'Need to select a superdir first.',
                                    ['Need to select a superdir first.']))
        self.add_parameter(Param('Static length [s]', 'static_length', 3, 1, 10))
        self.add_parameter(Param('Marker length [rel %]', 'marker_perc_length',\
                10., 0., 100.))

        self.add_parameter(Param('Number of Depths', 'num_depths', 3, 1, 21))
        self.add_parameter(Param('rel. Shift of Marker []', 'marker_shift_frac',\
                0.3, 0., 1.))

        self.add_parameter(Param('Number of Time Shifts', 'num_time_shifts',\
                6, 0, 30))

        self.add_parameter(Param('rel. Shift [%]', 'perc_of_shift', 20., 0., 100.))
        self.add_parameter(Param('Rise Time [s]', 'rise_time', 1., 0.5, 5.))
        self.add_trigger('Load Misfit Setup', self.load_misfit_setup) 
        self.add_trigger('Load Default Setup', self.load_setup) 
        self.add_trigger('Generate Markers', self.generate_markers) 
        self.add_trigger('Select Store Directory', self.setup_id_choice)
        self.add_trigger('Save', self.save)
        self.set_live_update(False)
        self.test_case_setup = None
        self.phase_ids_start = 'p|P'

    def call(self):
        '''Main work routine of the snuffling.'''

        self.cleanup()
        if not self.test_case_setup:
            self.fail(fail_message)
        self.active_event, self.stations = self.get_active_event_and_stations()

        #depths = self.test_case_setup.depths
        self.targets = du.stations2targets(self.stations)
        self.reference_source = du.event2source(active_event, 'DC')
        sources = du.test_event_generator(self.reference_source, depths)

        stf = [[0., self.rise_time],[0.,1.]]
        
        # TODO: Qt4 wie directory waehlen fuer engine dirs
        engine = LocalEngine(store_superdirs=store_dirs)
        test_case_setup = TestCaseSetup(reference_source=self.reference_source,
                                        sources=sources,
                                        targets=self.targets,
                                        engine=self.engine,
                                        store_id=self.store_id,
                                        misfit_setup=self.misfit_setup,
                                        source_time_function=sft,
                                        number_of_time_shifts=self.num_time_shifts,
                                        percentage_of_shift=self.perc_of_shift,
                                        phase_ids_start=self.phase_ids_start,
                                        static_length=self.static_length,
                                        marker_perc_length=self.marker_perc_length,
                                        marker_shift_frac=self.marker_shift_frac,
                                        depths=num_depths)


        test_case = core.TestCase(self.test_case_setup)
        core.Doer(test_case)

        tmp_out_dir = self.tempdir()

        TestCase.yaml_dump(pjoin(tmp_out_dir, 'derec_results.yaml'))
        TestCase.yaml_dump_setup(pjoin(tmp_out_dir, 'derec_setup.yaml'))

    def setup_id_choice(self):
        store_ids = self.input_filename(caption='Select a store') 
        self.store_id_choices.set_choices(store_ids)

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
                self.test_case_setup.number_of_time_shifts)
        self.set_parameter('num_depths', \
                len(self.test_case_setup.depths))
        self.set_parameter('rise_time', \
                self.test_case_setup.source_time_function[0][1])
        self.set_parameter('marker_perc_length', \
                self.test_case_setup.marker_perc_length)
        self.set_parameter('marker_shift_frac', \
                self.test_case_setup.marker_shift_frac)
        self.phase_ids_start = self.test_case_setup.phase_ids_start

    def save(self):
        self.output_filename('Save Results')

    def generate_markers(self):
        #if not self.test_case_setup:
        #    self.fail(fail_message)

        self.active_event, self.stations = self.get_active_event_and_stations()
        self.targets = du.stations2targets(self.stations)
        print self.active_event
        self.reference_source = du.event2source(self.active_event, 'DC')
        #fn = pjoin(derec_home, 'derec', 'test_case_setup.yaml')
        #f = open(fn,'r')
        #print 'loadl'
        #print f
        #self.test_case_setup = load_string(f.read())
        #print 'done'
        #f.close()
        #print self.test_case_setup


        # need to improve the horrible double line
        #store = self.test_case_setup.store_id 
        engine.regularize()
        print 'got it'
        markers = du.chop_ranges(self.reference_source,
                       self.targets,
                       engine.get_store(store_id),
                       self.phase_ids_start,
                       perc=self.marker_perc_length,
                       t_shift_frac=self.marker_shift_frac,
                       use_cake=True)

        print markers
        self.add_markers(markers.values()[0].values())

                
def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ Derec() ]


