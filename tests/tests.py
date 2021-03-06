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

    tpile_1 = list(pile.make_pile(pjoin(mseeds, 'testFiles/'), show_progress=False).iter_all())
    tpile_2 = list(pile.make_pile(pjoin(mseeds, 'testFiles_2/'), show_progress=False).iter_all())

    random_pile_test = list(pile.make_pile(pjoin(mseeds, 'testRandomTestFiles/'), show_progress=False).iter_all())
    random_pile_reference = list(pile.make_pile(pjoin(mseeds, 'referenceRandomTestFiles/'), show_progress=False).iter_all())

    t_min = util.str_to_time('2010-01-01 22:00:00')
    data1 = np.array([0, 0, 0, 0, 1, 0, 0, 0])
    data2 = np.array([0, 0, 0, 0, -1, 0, 0, 0])
    ttrace_1 = trace.Trace(network='1', station='TEST', channel='z', deltat=0.5, ydata=data1)
    ttrace_2 = trace.Trace(network='2', station='TEST', channel='z', deltat=0.5, ydata=data2)

    def test_misfit_by_samples(self):
        data = np.random.random(100)
        ttrace = trace.Trace(network='1', station='T1', channel='z', deltat=0.5, ydata=data)

        equal_traces_list = [ttrace.get_ydata(), ttrace.get_ydata()]
        self.assertEqual(misfit_by_samples(equal_traces_list), 0,
                         'misfit of same traces is zero')

    def test_misfit_by_samples_squared(self):
        equal_neg_traces_list = [self.ttrace_1.get_ydata(), self.ttrace_2.get_ydata()]
        self.assertEqual(misfit_by_samples(equal_neg_traces_list, square=True), 0,
                         'misfit of squared traces is zero if one trace is negative of the other')

    def test_find_matching_traces(self):
        self.assertEqual(len(find_matching_traces([list(self.tpile_1)], test_list=self.tpile_2)), 12,
                         'did not find all 12 of 12 pairs of traces')

    def test_time_domain_misfit_equal_piles(self):
        self.assertEqual(time_domain_misfit(reference_pile=[self.tpile_1],
                                            test_list=self.tpile_2), 0,
                         'Time Domain Misfit of equal piles is not 0!')

    def test_frequency_domain_misfit_equal_piles(self):
        self.assertEqual(frequency_domain_misfit(reference_pile=[self.tpile_1],
                                                 test_list=self.tpile_2), 0,
                         'Frequency Domain Misfit of equal test_piles is NOT 0')

    def test_time_domain_misfit_UNequal_piles(self):
        self.assertNotEqual(time_domain_misfit(reference_pile=[self.random_pile_reference],
                                               test_list=self.random_pile_test), 0,
                            'Time Domain Misfit of UNequal traces is not 0')

    def test_frequency_domain_misfit_UNequal_piles(self):
        self.assertNotEqual(frequency_domain_misfit(reference_pile=[self.random_pile_reference],
                                                    test_list=self.random_pile_test), 0,
                            'Freq Dom Misfit of unequal piles is 0 but should not be 0')

    def test_equalize_sampling_rate(self):
        '''

        :return:
        '''
        data1 = np.random.random(100)
        data2 = np.random.random(100)
        ttrace_1 = trace.Trace(network='1', station='T1', channel='z', deltat=0.5, ydata=data1)
        ttrace_2 = trace.Trace(network='1', station='T2', channel='z', deltat=0.01, ydata=data2)
        downsample_if_needed([[ttrace_1, ttrace_2]])
        self.assertEqual(ttrace_1.deltat, ttrace_2.deltat, 'Equalization of sampling rates unsuccessful!')
        self.assertEqual(ttrace_1.deltat, 0.6, 'new sampling rate of ttrace_1 wrong')
        self.assertEqual(ttrace_2.deltat, 0.6, 'new sampling rate of ttrace_2 wrong')

    def test_chop_to_same_sample_length(self):
        data1 = np.random.random(100)
        data2 = np.random.random(80)
        [data1, data2] = chop_longer_samples([data1, data2])
        self.assertEqual(np.shape(data1), np.shape(data2), 'shape after chopping not equal')

#
#class TestPyrocko(unittest.TestCase):
#    def setUp(self):
#        self.tpile_1 = pile.make_pile('../mseeds/testFiles', show_progress=False)
#        self.tpile_2 = pile.make_pile('../mseeds/testFiles', show_progress=False)
#
#
#    def test_data_in_both_piles(self):
#        self.assertTrue(noDatalessTracesInPile(self.tpile_1),
#                        'some dataless traces in test_pile_1')
#        self.assertTrue(noDatalessTracesInPile(self.tpile_2),
#                        'some dataless traces in test_pile_2')
#
#    def test_add_trace_to_pile(self):
#        ydata = np.array([1,2,3,4,5,51,2,3,2,2])
#        tt = trace.Trace(station='test', ydata=ydata)
#
#        tp = pile.Pile()
#        tp.add(tt)
#        self.assertTrue(isinstance(tp, pile.Pile))
#        self.assertTrue(isinstance(tt, trace.Trace))

if __name__ == '__main__':
    unittest.main()


