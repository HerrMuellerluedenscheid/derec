from pyrocko.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko import trace, gf
from derec import derec_utils as du
from derec import core
from derec import yaml_derec
from pyrocko.guts import load_string, String, Float, Int, Dict
from pyrocko.guts_array import Array
import numpy as num
import tempfile
import os

pjoin = os.path.join
fail_message = 'Need to load a setup, first'

class Derec(Snuffling):
    
    '''
    '''

    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Depth Relocation')
        self.add_parameter(Choice('Store: ', 'store_id', 
                                    'Need to select a superdir first.',
                                    'Need to select a superdir first.'))
        self.add_parameter(Param('Number of Depths', 'depths', 3, 1, 21))
        self.add_parameter(Param('Number of Time Shifts', 'num_time_shifts',\
                9, 1, 21))
        self.add_trigger('Load Default Setup', self.load_setup) 
        self.add_trigger('Generate Markers', self.generate_markers) 
        self.add_trigger('Select Store Directory', self.setup_id_choice)
        self.add_trigger('Save', self.save)
        self.set_live_update(False)
        self.test_case_setup = None

    def call(self):
        '''Main work routine of the snuffling.'''

        self.cleanup()
        if not self.test_case_setup:
            self.fail(fail_message)

        self.active_event, self.stations = self.get_active_event_and_stations()

        depths = self.test_case_setup.depths
        self.targets = du.stations2targets(self.stations)
        self.reference_source = du.event2source(active_event, 'DC')
        sources = du.test_event_generator(reference_source, depths)

        test_case = core.TestCase(self.test_case_setup)
        core.Doer(test_case)

        tmp_out_dir = self.tempdir()

        TestCase.yaml_dump(pjoin(tmp_out_dir, 'derec_results.yaml'))
        TestCase.yaml_dump_setup(pjoin(tmp_out_dir, 'derec_setup.yaml'))

    def setup_id_choice(self):
        store_ids = self.input_filename(caption='Select a store') 
        self.store_id_choices.set_choices(store_ids)

    def load_setup(self):
        fn = self.input_filename(caption='Select a setup')
        f = open(fn,'r')
        self.test_case_setup = load_string(f.read())
        f.close()
        
        self.set_parameter('num_time_shifts', \
                self.test_case_setup.number_of_time_shifts)

    def save(self):
        self.output_filename('Save Results')

    def generate_markers(self):
        if not self.test_case_setup:
            self.fail(fail_message)

        self.active_event, self.stations = self.get_active_event_and_stations()
        self.targets = du.stations2targets(self.stations)
        self.reference_source = du.event2source(self.active_event, 'DC')

        # need to improve the horrible double line
        self.markers = du.chop_ranges(self.reference_source,
                       self.targets,
                       self.test_case_setup.engine.get_store(
                           self.test_case_setup.store_id),
                       self.test_case_setup.phase_ids_start,
                       perc=self.test_case_setup.marker_perc_length,
                       t_shift_frac=self.test_case_setup.marker_shift_frac,
                       use_cake=True)

        self.add_markers(self.markers.values()[0].values())

                
def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ Derec() ]


