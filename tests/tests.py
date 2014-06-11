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
        num_sources = 100
        sources_lists = du.make_lots_of_test_events(self.source, [1000,2000,3000], 
                {'strike':10,
                'dip':10,
                'rake':10}, 
                num_sources)

        assert len(sources_lists)==num_sources
        for i in range(num_sources-1):
            for n in range(len(sources_lists[0])):
                assert sources_lists[0][n].dip!=sources_lists[i+1][n].dip

    def test_random_event_generation_range(self):
        num_sources = 100
        sources_lists = du.make_lots_of_test_events(self.source, [1000,2000,3000], 
                {'strike':10,
                'dip':10,
                'rake':[10,20]}, 
                num_sources)
        import pdb
        pdb.set_trace()
        assert len(sources_lists)==num_sources
        for i in range(num_sources-1):
            for n in range(len(sources_lists[0])):
                assert sources_lists[0][n].dip!=sources_lists[i+1][n].dip

        

if __name__ == '__main__':
    unittest.main()


