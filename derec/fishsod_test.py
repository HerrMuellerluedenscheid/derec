import unittest
import numpy as np
from pyrocko import pile, trace, util
from fishsod_utils import *


def noDatalessTracesInPile(p):
    for t in p.iter_traces(load_data=True):
        try:
            t.get_ydata()
        except trace.NoData:
            return False
    return True


class TestFSU(unittest.TestCase):

    tpile_1 = pile.make_pile('testFiles/', show_progress=False)
    tpile_2 = pile.make_pile('testFiles_2/', show_progress=False)

    random_pile_test = pile.make_pile('testRandomTestFiles/', show_progress=False)
    random_pile_reference = pile.make_pile('referenceRandomTestFiles/', show_progress=False)

    traces = []
    for tr in tpile_1.iter_traces(load_data=True):
        traces.append(tr.copy(data=True))

    t_min = util.str_to_time('2010-01-01 22:00:00')
    data1 = np.array([0, 0, 0, 0, 1, 0, 0, 0])
    data2 = np.array([0, 0, 0, 0, -1, 0, 0, 0])
    ttrace_1 = trace.Trace(network='1', station='TEST', channel='z', deltat=0.5, ydata=data1)
    ttrace_2 = trace.Trace(network='2', station='TEST', channel='z', deltat=0.5, ydata=data2)

    def test_make_pile(self):
        test_pile_4 = pile.make_pile('testFiles/', show_progress=False)
        self.assertNotEqual(test_pile_4, self.tpile_1,
                            'Test Piles from different directories are equal.')
        self.assertNotEqual(test_pile_4.all()[0].get_ydata(), None,
                            'Data of first trace is None.')

    def test_misfit_by_samples(self):
        equal_traces_list = [self.traces[1].get_ydata(), self.traces[1].get_ydata()]
        self.assertEqual(misfit_by_samples( equal_traces_list ), 0,
                         'misfit of same traces is zero')

    def test_misfit_by_samples_squared(self):
        eqal_neg_traces_list = [self.ttrace_1.get_ydata(), self.ttrace_2.get_ydata()]
        self.assertEqual(misfit_by_samples( eqal_neg_traces_list, square=True), 0,
                         'misfit of squared traces is zero if one trace is negative of the other')

    def test_find_matching_traces(self):
        self.assertEqual(len(find_matching_traces(self.tpile_1, self.tpile_2)), 12,
                         'did not find all 12 of 12 pairs of traces')

    def test_time_domain_misfit_equal_piles(self):
        self.assertTrue(noDatalessTracesInPile(self.tpile_1),
                        'Dataless traces in pile found')
        self.assertTrue(noDatalessTracesInPile(self.tpile_2),
                        'Dataless traces in pile found')
        self.assertEqual(time_domain_misfit(self.tpile_1, self.tpile_2), 0,
                         'Time Domain Misfit of equal piles is not 0!')

    def test_frequency_domain_misfit_equal_piles(self):
        # evtl muss da getapert werden.
        self.assertTrue(noDatalessTracesInPile(self.tpile_1),
                        'Dataless traces in pile found')
        self.assertTrue(noDatalessTracesInPile(self.tpile_2),
                        'Dataless traces in pile found')
        self.assertEqual(frequency_domain_misfit(self.tpile_1, self.tpile_2), 0,
                         'Frequency Domain Misfit of equal test_piles is NOT 0')

    def test_time_domain_misfit_of_random_traces(self):
        random_pile_test = pile.make_pile('mseeds/testRandomTestFiles/', show_progress=False)
        random_pile2_reference = pile.make_pile('mseeds/referenceRandomTestFiles/', show_progress=False)

        random_list_test=[tr for tr in random_pile_test.iter_traces(load_data=True)]
        self.assertNotEqual(time_domain_misfit(reference_pile=random_pile2_reference,
                                               test_list=random_list_test), 0,
                            'MF of 2piles with random traces is 0, should not be 0')

    def test_time_domain_misfit_UNequal_piles(self):
        self.assertNotEqual(time_domain_misfit(self.random_pile_reference, self.random_pile_test), 0,
                            'Time Domain Misfit of UNequal traces is not 0')

    def test_frequency_domain_misfit_UNequal_piles(self):
        # evtl muss da getapert werden.
        self.assertNotEqual(frequency_domain_misfit(self.random_pile_reference, self.random_pile_test), 0,
                            'Freq Dom Misfit of unequal piles is 0 but should not be 0')


class TestPyrocko(unittest.TestCase):
    def setUp(self):
        self.tpile_1 = pile.make_pile('testFiles', show_progress=False)
        self.tpile_2 = pile.make_pile('testFiles', show_progress=False)


    def test_data_in_both_piles(self):
        self.assertTrue(noDatalessTracesInPile(self.tpile_1),
                        'some dataless traces in test_pile_1')
        self.assertTrue(noDatalessTracesInPile(self.tpile_2),
                        'some dataless traces in test_pile_2')

    def test_add_trace_to_pile(self):
        ydata = np.array([1,2,3,4,5,51,2,3,2,2])
        tt = trace.Trace(station='test', ydata=ydata)

        tp = pile.Pile()
        tp.add(tt)
        self.assertTrue(isinstance(tp, pile.Pile))
        self.assertTrue(isinstance(tt, trace.Trace))

if __name__ == '__main__':
    unittest.main()


