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
                {('strike', 'dip', 'rake'):10}, 
                num_sources)

        assert len(sources_lists)==num_sources
        for i in range(num_sources-1):
            for n in range(len(sources_lists[0])):
                assert sources_lists[0][n].dip!=sources_lists[i+1][n].dip

    def test_random_event_generation_range(self):
        num_sources = 100
        source = du.clone(self.source)
        source.strike = 130
        source.dip = 88
        source.rake = -119
        sources_lists = du.make_lots_of_test_events(source, [1000,2000,3000], 
                {('strike', 'dip', 'rake'):30}, 
                num_sources,
                func='normal')
        assert len(sources_lists)==num_sources
        for i in range(num_sources-1):
            for n in range(len(sources_lists[0])):
                assert sources_lists[0][n].dip!=sources_lists[i+1][n].dip

        dips =[]
        rakes = []
        strikes = []
        angles = []
        for sl in sources_lists:
            s = sl[0]
            strikes.append(s.strike)
            rakes.append(s.rake)
            dips.append(s.dip)
            a = s.pyrocko_moment_tensor().angle(source.pyrocko_moment_tensor())
            angles.append(a)
        import matplotlib.pyplot as plt
        f, axs = plt.subplots(5)
        axs[0].hist(strikes, 25)
        axs[1].hist(dips, 25)
        axs[2].hist(rakes, 25)
        axs[3].hist(angles, 25)

        gauss = [num.random.normal(0,10) for i in xrange(num_sources)]
        axs[4].hist(num.abs(gauss), 25)

        #plt.show()

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
        Mfinal, return_scaling = du.L2_norm(can_tr, ref_tr)
        assert Mfinal.values()[0]==0.

        # add random noise:
        for s,t,tr in du.iter_dict(ref_tr):
            du.add_random_noise_to_trace(tr,A=max(abs(tr.get_ydata()))/10)
        Mfinal, return_scaling = du.L2_norm(can_tr, ref_tr)
        assert Mfinal.values()[0]!=0.
        
        for s,t,tr in du.iter_dict(ref_tr):
            tr.ydata=num.ones(len(tr.ydata))
        for s,t,tr in du.iter_dict(can_tr):
            tr.ydata=num.zeros(len(tr.ydata))
        
        Mfinal, return_scaling = du.L2_norm(can_tr, ref_tr)

        # use scaling c
        for s,t,tr in du.iter_dict(ref_tr):
            tr.ydata=num.ones(len(tr.ydata))
        for s,t,tr in du.iter_dict(can_tr):
            tr.ydata=num.ones(len(tr.ydata))
            tr.ydata*=0.5

        c = 2.
        Mfinal, return_scaling = du.L2_norm(can_tr, ref_tr, scaling=[c])
        assert Mfinal.values()[0]==0.

    def test_snr_hisignal(self):
        for i in range(10):
            no = num.random.uniform(-1, 1, 100000)
            si = num.random.uniform(-10, 10, 100000)
            ts = trace.Trace(ydata=si)
            tn = trace.Trace(ydata=no)
            assert (core.snr(ts, tn)-20)<0.1

    def test_snr_hinoise(self):
        """ difference should be 50 db according to 
        http://en.wikipedia.org/wiki/20_log_rule"""
        for i in range(10):
            no = num.random.normal(0, 1, 100000)
            si = num.random.normal(0, 0.003162, 100000)
            ts = trace.Trace(ydata=si)
            tn = trace.Trace(ydata=no)
            assert (core.snr(ts, tn)+50.)<0.1

    def test_snr_same_level(self):
        for i in range(10):
            no = num.random.normal(0, 1, 100000)
            si = num.random.normal(0, 1, 100000)
            ts = trace.Trace(ydata=si)
            tn = trace.Trace(ydata=no)
            assert num.abs(core.snr(ts, tn))<0.1


        

if __name__ == '__main__':
    unittest.main()


