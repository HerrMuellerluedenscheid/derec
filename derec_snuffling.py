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
        self.add_trigger('Select Store Directory', self.setup_id_choice)
        self.add_trigger('Save', self.save)
        self.set_live_update(False)
        self.test_case_setup = None

    def call(self):
        '''Main work routine of the snuffling.'''
        if not self.test_case_setup:
            self.warn('Need to select a setup, first.')

        self.cleanup()

        # als window:
        #store_superdir = ''
        # und dann als drop down menu
        # note: jedes Feld sollte als yaml Feld gebaut werden um snuffler zustand
        #    fuers naechset Mal zu speichern.
        
        # kann auch als drop down, also local oder remote engine
        #engine = LocalEngine(store_superdirs=store_superdirs)
        #store_ids = engine.get_store_ids()

        event, stations = self.get_active_event_and_stations()
        self.test_case_setup.targets = du.stations2targets(stations)
        print event
        test.test_case_setup.reference_source = du.event2source(event, 'DC')
        
        #n_z = 5
        #z_offset = 1000
        #depths = num.linspace(source.depth-z_offset,
        #                    source.depth+z_offset, 
        #                    n_z)

        test_case = core.TestCase(test_case_setup)
        core.Doer(test_case)

        print test_case.misfits
        tmp_out_dir = self.tempdir()

        TestCase.yaml_dump(pjoin(tmp_out_dir, 'derec_results.yaml'))
        TestCase.yaml_dump_setup(pjoin(tmp_out_dir, 'derec_setup.yaml'))

    def setup_id_choice(self):
        self.store_ids = self.input_filename(caption='Select a store') 
        self.store_id_choices.set_choices(self.store_ids)

    def load_setup(self):
        fn = self.input_filename(caption='Select a setup')
        f = open(fn,'r')
        self.test_case_setup = load_string(f.read())
        f.close()
        self.set_parameter('num_time_shifts', \
                self.test_case_setup.number_of_time_shifts)

    def save(self):
        self.output_filename('Save Results')
                
def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ Derec() ]


