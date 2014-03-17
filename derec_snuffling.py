from pyrocko.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko import trace, gf
from derec import derec_utils as du
from derec import Core
import numpy as num

class Derec(Snuffling):
    
    '''
    '''

    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Depth Relocation')
        self.store_id_choice = self.Choice('Store: ', 'store_id', 
                                    'Need to select a superdir first.')
        self.add_trigger('Select Store Directory', setup_id_choice)
        self.set_live_update(False)

    def call(self):
        '''Main work routine of the snuffling.'''
        
        self.cleanup()
        viewer = self.get_viewer()
        pile = self.get_pile()
        
        # als window:
        store_superdir = ''
        # und dann als drop down menu
        # note: jedes Feld sollte als yaml Feld gebaut werden um snuffler zustand
        #    fuers naechset Mal zu speichern.
        
        # kann auch als drop down, also local oder remote engine
        engine = LocalEngine(store_superdirs=store_superdirs)
        store_ids = engine.get_store_ids()

        markers = viewer.selected_markers()
        event, stations = self.get_active_event_and_stations()
        
        n_z = 5
        z_offset = 100 
        depths = num.arange(source.depth-z_offset,
                            source.depth+z_offset, 
                            n_z)

        core = Core(markers=markers, 
                    stations=stations, 
                    misfit_setup=misfit_setup, 
                    test_case_setup=test_case_setup)

    def setup_id_choice(self):
        self.store_ids = self.input_filename(caption='Select a store') 
        self.store_id_choices.set_choices(self.store_ids)
                
def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ Derec() ]


