import unittest
from os import path
from pyrocko import trace, pile, util
from pyrocko.gf.seismosizer import *
import derec.derec_utils as du
import numpy as np
import os

pjoin = path.join
dhome = path.dirname(__file__)
derec_home = os.environ['DEREC_HOME']
mseeds = path.abspath(pjoin(dhome, '../mseeds'))


def no_dataless_traces_in_pile(p):
    for t in p.iter_traces(load_data=True):
        try:
            t.get_ydata()
        except trace.NoData:
            return False
    return True


class TestFSU(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestFSU, self).__init__(*args, **kwargs)
        self.store_id = 'doctar_mainland_20Hz'
        self.target = Target()
        self.source = DCSource()
        
    
    def test_random_event_generation(self):
        sources_lists = du.make_lots_of_test_events(self.source, [1000,2000,3000], 
                {'strike':10,
                'dip':10,
                'rake':10}, 
                10000)

        assert len(sources_lists)==10000
        for i in range(10000-1):
            for n in range(len(sources_lists[0])):
                assert sources_lists[0][n].dip!=sources_lists[i+1][n].dip

    def test_generate_markers(self):
        store_id = 'castor'
        target = Target(lat=10., lon=10., store_id=store_id)
        source = DCSource(lat=10., lon=10.1, depth=5000)
        engine = LocalEngine(store_superdirs=[derec_home+'/fomostos'])
        #self.active_event, self.stations = self.get_active_event_and_stations()
 
        #if not self.targets:
        #    self.targets = du.stations2targets(self.stations, \
        #            self.store_id_choice)
 
        
        #    self.reference_source = DCSource.from_pyrocko_event(self.active_event)
        #    self.__reference_source_dict = self.reference_source.__dict__
 
        tmins = [] 
        #for sl in [2.,3.]:
        a, c = du.chop_ranges(source,
                               [target],
                               engine.get_store(store_id),
                               ['p', 'P'],
                               return_cache=True,
                               cache=True,
                               perc=0.0001,
                               static_length=1,
                               t_shift_frac=0.5,
                               use_cake=True)
        
        del c #c.flush()
        self._ma = a.values()[0].values()
        tmins.append(self._ma[0].tmin)
        #assert tmins[0]!=tmins[1]

        print a.values()[0].values()[0].tmin
 
        b = du.chop_ranges(source,
                               [target],
                               engine.get_store(store_id),
                               ['p'],
                               return_cache=False,
                               cache=False,
                               perc=0.0001,
                               static_length=2,
                               t_shift_frac=0.5,
                               use_cake=True)
        m = b.values()[0].values()[0]
        print m.tmin
        print m.tmax-m.tmin

        

if __name__ == '__main__':
    unittest.main()


