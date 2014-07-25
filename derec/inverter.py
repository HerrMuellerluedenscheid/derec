from collections import defaultdict
from derec.yaml_derec import TestCaseSetup
from derec.core import TestCase
from derec import derec_utils as du
from derec import optics
from pyrocko.gf import meta
import forkmap
import time
import numpy as num
map = forkmap.map


def chop_ranges_reformatted(chops):
    chops_reform = defaultdict()
    for s,t_m in chops.items():
        for t, m in t_m.items():
            d = s.distance_to(t)
            z = s.depth
            chops_reform[(d,z)] = m

    return chops_reform


class Inverter(object):
    def __init__(self, test_case, debug=False):
        self.debug = debug
        self.verbose = False
        self.test_case = test_case
        self.setup = test_case.test_case_setup
        self.setup_dict = self.setup.__dict__

        self.reference_seismograms = test_case.raw_references
        self.last_best_source = self.setup.reference_source
        self.last_best_mf = 999.
        self.ref_markers = self.test_case.reference_markers
        self.noise = None
        self.steps = 5
        self.tests_per_step = 10
        self.best_steps = []
        self.focus = None


    def get_new_setup(self):
        return TestCaseSetup(**self.setup_dict)

    def get_new_test_case(self):
        return TestCase(self.setup)

    def get_result_str(self):
        'Inversion results: \n'
        outstr = ''
        for mf,s in self.best_steps:
            outstr += 'M=%s at strike = %s,  dip = %s,  rake = %s \n'%(mf, s.strike, s.dip, s.rake)
        return outstr
    
    def get_results_unformatted(self):
        outstr = ''
        for mf,s in self.best_steps:
            outstr += '%s %s %s %s \n'%(mf, s.strike, s.dip, s.rake)
        return outstr


    def print_results(self):
        print self.get_result_str()
    
    def setup_cases(self, location_test_sources):
        outp = []
        for s in location_test_sources:
            test_case_setup = self.get_new_setup()
            test_case_setup.sources = s
            test_case = TestCase(test_case_setup)
            test_case.set_raw_references(self.reference_seismograms)
            test_case.set_reference_markers(self.ref_markers)
            outp.append(test_case)

        return outp



class SimulAn(Inverter):
    def __init__(self, test_case, debug=False):
        super(SimulAn, self).__init__(test_case, debug)

    def accept(self, bestmf, lastbestmf, num_inversions, i_total, circles_before_break,
            debug=False):
        if bestmf<lastbestmf:
            return 1 
        e = 3.

        T = num.linspace(0,1, circles_before_break)
        T = N*num.exp(-e*T)

        T_c = num.random.random(1)
        if T_c > T[i_total]:
            return 0 
        else:
            return 1 


class FocusMonteCarlo(Inverter):
    def __init__(self, *arg, **kwargs):
        super(FocusMonteCarlo, self).__init__(*arg, **kwargs)
        pass


    def set_focus(self, focus=None):
        if focus:
            self.focus = focus
            self.steps = len(focus)-1
        elif self.focus==None:
            foci = 90*num.exp(-3.5*num.linspace(0,1, self.steps))
            self.focus = num.ones(self.steps+1)*360
            self.focus[1::] = foci

    def run(self, nworkers=1, verbose=False):
        self.verbose = verbose
        for i in range(self.steps+1):
            location_test_sources_lists = self.get_test_sources_lists(i)
            cases = self.setup_cases(location_test_sources_lists)
            map(self.run_step, cases, n=nworkers)
            time.sleep(0.1)
            self.evaluate_last(cases)

    def evaluate_last(self, cases):
        minmf = self.last_best_mf
        for c in cases:
            s, m = c.best_source_misfit()
            if m<minmf:
                minmf=m
                bests=s
            
        if minmf<self.last_best_mf:
            self.last_best_mf = minmf
            self.last_best_source = bests

        self.best_steps.append((self.last_best_mf, self.last_best_source))

    def get_test_sources_lists(self, step):
        self.set_focus()
        focus = self.focus[step]
        if self.debug or self.verbose: print 'new focus: ', focus
    
        tevents = du.make_lots_of_test_events(self.last_best_source,
                                           self.setup.depths,
                                           {('strike', 'dip', 'rake'):focus},
                                           self.tests_per_step,
                                           isfirst=(step==0),
                                           func='uniform')

        if self.debug:
            print 'checking strike dip rake setup...'
            optics.check_locations(tevents)
        return tevents

    def run_step(self, test_case):
        try:
            test_case.process(verbose=False, use_cake=True)
        except meta.OutOfBounds:
            pass







    
