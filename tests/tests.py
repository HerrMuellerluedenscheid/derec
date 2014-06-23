import unittest
from os import path
from pyrocko import trace, pile, util
from pyrocko.gf.seismosizer import *
import derec.derec_utils as du
import numpy as num
import os
import derec.core as core

pjoin = path.join
derec_home = os.environ['DEREC_HOME']
mseeds = path.abspath(pjoin(derec_home, 'mseeds'))
fomostos = path.abspath(pjoin(derec_home, 'fomostos'))


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
        self.engine = LocalEngine(store_superdirs=[fomostos])
        
    
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
        assert len(sources_lists)==num_sources
        for i in range(num_sources-1):
            for n in range(len(sources_lists[0])):
                assert sources_lists[0][n].dip!=sources_lists[i+1][n].dip

    def test_my_L2Norm(self):
        lon = 10.
        lat = 10.
        channels = 'NEZ'

        targets = [Target(lon=lon,
                          lat=lat, 
                          codes=('Y7', 'L001', '', c), 
                          store_id=self.store_id) for c in channels]
        
        source = DCSource(lon=lon, lat=lat+0.2, depth=5000)
        ref_tr = core.make_reference_trace(source, targets, self.engine)
        can_tr = core.make_reference_trace(source, targets, self.engine)

        assert du.L2_norm(can_tr, ref_tr).values()[0]==0.

        # add random noise:
        for s,t,tr in du.iter_dict(ref_tr):
            du.add_random_noise_to_trace(tr,A=max(abs(tr.get_ydata()))/10)
        assert du.L2_norm(can_tr, ref_tr).values()[0]!=0 
        
        for s,t,tr in du.iter_dict(ref_tr):
            tr.ydata=num.ones(len(tr.ydata))
        for s,t,tr in du.iter_dict(can_tr):
            tr.ydata=num.zeros(len(tr.ydata))
        
        assert du.L2_norm(can_tr, ref_tr).values()[0]==1

        # use scaling c
        for s,t,tr in du.iter_dict(ref_tr):
            tr.ydata=num.ones(len(tr.ydata))
        for s,t,tr in du.iter_dict(can_tr):
            tr.ydata=num.ones(len(tr.ydata))
            tr.ydata*=0.5

        c = 2.
        assert du.L2_norm(can_tr, ref_tr, scaling=c).values()[0]==0

        

if __name__ == '__main__':
    unittest.main()


