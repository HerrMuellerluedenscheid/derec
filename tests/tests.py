import unittest
from os import path
from pyrocko import trace, pile, util
from derec.derec_utils import *
import numpy as np

pjoin = path.join
dhome = path.dirname(__file__)
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
        self.store_id = 'castor'
        self.target = Target()
        self.source = DCSource()
    
    def test_random_event_generation(self):
        sources_lists = make_lots_of_test_events(self.source, [1000,2000,3000], 
                {'strike':10,
                'dip':10,
                'rake':10}, 
                10000)

        assert len(sources_lists)==10000
        for i in range(10000-1):
            for n in range(len(sources_lists[0])):
                assert sources_lists[0][n].dip!=sources_lists[i+1][n].dip


        

if __name__ == '__main__':
    unittest.main()


