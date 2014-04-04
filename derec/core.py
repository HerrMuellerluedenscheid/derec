from yaml_derec import *
from pyrocko.gf import *
from pyrocko import model, gui_util, trace, moment_tensor, io
from collections import defaultdict
from scipy.signal import butter
from scipy.ndimage import zoom
from pyrocko.guts import Object, Float, Int, String, Complex, Tuple, List, load_string, Dict
from pyrocko.guts_array import Array

import time
import matplotlib.lines as pltlines
import progressbar
import os
import derec_utils as du
import numpy as num
import copy

pjoin = os.path.join
km = 1000.


def equal_attributes(o1, o2):
    '''
    Return true if two objects are equal as for their attributes. 
    '''
    return o1.__dict__ == o2.__dict__


def set_refine_parameter(ref_event, **kwargs):
    '''
    Returns dict. Key is the value of **kwargs. 
    '''
    events = {}
    for k, vals in kwargs.iteritems():
        for val in vals:
            # mit event.copy ersetzen?
            event_copy = copy.copy(ref_event)
            exec('event_copy.%s=%s' % (k, val))
            events[val]=event_copy

    return events


def make_reference_trace(source, targets, engine):
    if not isinstance(source, list):
        source = [source]

    response = engine.process(
            sources=source,
            targets=targets)
    return response
    

class Doer():
    def __init__(self, test_case):

        test_case.request_data()

        #print 'source location: ', test_case.ref_source
        print('test data marker....')
        extended_test_marker = du.chop_ranges(test_case.sources,
                                              test_case.targets,
                                              test_case.store,
                                              test_case.test_case_setup.phase_ids_start,
                                              perc=1.0,
                                              t_shift_frac=0.3)
        
        test_case.set_candidates_markers( extended_test_marker )

        print('chopping ref....')
        test_case.references = du.chop_using_markers(
                                test_case.raw_references, 
                                test_case.reference_markers, 
                                inplace=False)

        test_case.apply_stf(test_case.test_case_setup.source_time_function)

        print('chopping cand....')
        test_case.candidates = du.chop_using_markers(
                                test_case.raw_candidates, 
                                extended_test_marker, 
                                inplace=False)

        du.calculate_misfit(test_case)



class TestCase(Object):
    '''
    In one test case, up to 3 parameters can be modified
    '''
    def __init__(self, test_case_setup):
        self.test_case_setup = test_case_setup

        self.targets = test_case_setup.targets
        self.sources = test_case_setup.sources
        self.engine = test_case_setup.engine
        self.test_parameters = test_case_setup.test_parameters 
        self.store_id = test_case_setup.store_id

        self.raw_references = None
        self.processed_references = defaultdict(dict)
        self.references = {} 

        self.raw_candidates = None
        self.processed_candidates = defaultdict(dict)
        self.candidates= {}
        self.misfits = defaultdict
        self.misfit_setup = test_case_setup.misfit_setup
        self.channel_map = test_case_setup.channel_map    
            
    def request_data(self):
        print 'requesting data....'
        self.response = self.engine.process(status_callback=self.update_progressbar, 
                                sources=self.sources,
                                targets=self.targets)
        print 'finished'
        self.set_raw_candidates(du.response_to_dict(self.response))
    
    def set_raw_references(self, references):
        """
        references is a dictionary containing the reference seismograms
        """
        self.raw_references = references 

    def set_raw_candidates(self, candidates):
        """
        candidates is a dictionary containing the candidates seismograms
        """
        self.raw_candidates = candidates

    def set_reference_markers(self, markers):
        """
        Reference markers are the ones used as stencil to chop the reference
        traces.
        """
        self.reference_markers = markers

    def set_candidates_markers(self, markers):
        """
        candidates markers are the ones used as stencil to chop the candidates 
        traces.
        """
        self.candidates_markers = markers

    def extend_markers(self, markers, c):
        """
        Scale a markers extension with factor *c*
        """
        markers = du.extend_markers(markers, scaling_factor=c)
        return markers

    def get_lowest_misfit_data(self):
        """
        Return source, target and the value of the lowest misfit.
        """
        misfit_value = 999.
        for s, mf in self.misfits.iteritems():
            if mf < misfit_value:
                source = s
                misfit_value = mf
        
        return source, misfit_value

    @staticmethod
    def targets_nsl_of(targets=[]):
        """return a set of all network, station, location tuples contained
        in *targets*"""
        return set([t.codes[:3] for t in targets])

    def apply_stf(self, stf):
        """
        Apply source time function on candidates.
        """
        self.raw_candidates = du.apply_stf(self.raw_candidates, stf)

    @property
    def targets_nsl(self):
        """return a set of all network, station, location tuples contained
        in this Test Case"""
        return self.targets_nsl_of(self.targets)

    def set_misfit(self, misfits):
        self.misfits = misfits 

    def yaml_dump(self, fn=''):

        def convert_to_yaml_dict(_dict):
            outdict = defaultdict(dict)
            for source, target_o in _dict.iteritems():
                for target, o in target_o.iteritems():
                    if isinstance(o, trace.Trace):
                        yaml_o = yamlTrace(ydata=o.ydata,
                                tmin=o.tmin,
                                deltat=o.deltat, 
                                codes=o.nslc_id)
                    elif isinstance(o, gui_util.Marker):
                        yaml_o = yamlMarker(nslc_ids=o.nslc_ids,
                                tmin=o.tmin,
                                tmax=o.tmax,
                                kind=o.kind)

                    outdict[source][target] = yaml_o

            return outdict

        test_case_data = TestCaseData()
        test_case_data.references = convert_to_yaml_dict(self.references)
        test_case_data.candidates = convert_to_yaml_dict(self.candidates)
        test_case_data.processed_references = convert_to_yaml_dict(
                                                self.processed_references)
        test_case_data.processed_candidates = convert_to_yaml_dict(
                                                self.processed_candidates)

        test_case_data.test_case_setup = self.test_case_setup

        misfit_float_dict = dict(zip(self.misfits.keys(),
                                [float(i) for i in self.misfits.values()]))

        test_case_data.results = dict(misfit_float_dict)

        test_case_data.reference_markers = convert_to_yaml_dict(
                self.reference_markers)
        test_case_data.candidates_markers = convert_to_yaml_dict(
                self.candidates_markers)

        f = open(fn, 'w')
        f.write(test_case_data.dump())
        f.close()

    def make_t_shifts(self, trac, num_samples, perc):
        """
        :param trac: pyrocko.trace.Trace
        :param num_samples: number of time shifts
        :param perc: percentage of trace length to be shifted 
        :return: numpy array. 
        """
        t_shift_max = (trac.tmax - trac.tmin) / 100. * perc
        return num.linspace(-t_shift_max/2., t_shift_max/2, num_samples)

    def make_shifted_candidates(self, source, target):
        """
        returns shifted candidates.
        """
        shifted_candidates = []
        cand = self.candidates[source][target]
        t_shifts = self.make_t_shifts(cand,
                self.test_case_setup.number_of_time_shifts, 
                self.test_case_setup.percentage_of_shift)

        shifted_candidates = [cand.copy() for i in range(len(t_shifts))]
        map(lambda t,s: t.shift(s), shifted_candidates, t_shifts)
        return shifted_candidates

    @staticmethod
    def yaml_read(fn):
        '''
        Create a TestCase Object from file. 

        :param fn: (str) filename
        '''
        try:
            f = open(fn, 'r')
            tc = load_string(f.read())
        except:
            pass
        finally:
            f.close()
            raise
        
    def update_progressbar(self, a, b):
        try:
            self.progressbar.update(a)
        except AttributeError:
            self.progressbar = progressbar.ProgressBar(maxval=b).start()
            self.progressbar.update(a)

    def snuffle(self):
        trace.snuffle(self.seismograms)

    def numpy_it(self, **kwargs):
        '''
        '''
        if kwargs.get('order', False):
            self.xkey, self.ykey, self.zkey = kwargs['order']
        else:
            self.xkey, self.ykey, self.zkey = self.test_parameters.keys()

        self.num_array = num.empty(shape=(len(self.sources), 4))
        self.num_array[:] = num.NAN

        for i, s in enumerate(self.sources):
            self.num_array[i] = num.array([getattr(s, self.xkey),
                                           getattr(s, self.ykey),
                                           getattr(s, self.zkey),
                                           self.misfits[s]])

        self.num_array = self.num_array.T

    @staticmethod
    def get_sources_where(param_dict, sources=None):
        '''
        param_dict is something like {latitude:10, longitude:10}

        :returns: list of sources, matching the required parameters in the param_dict.
        '''
        return filter(lambda s: all(map(lambda x: abs(getattr(s, x[0])-x[1])<1e-7, \
                param_dict.iteritems())), sources)


    @staticmethod
    def lines_dict(traces_dict):
        """
        Create matplotlib.lines.Line2D objects from traces dicts.
        :return lines_dict: dict with lines
        """
        lines_dict = defaultdict(dict)
        for source, target, tr in TestCase.iter_dict(traces_dict):
            # this is ugly and shouldn't be necessary if get_xdata() worked:
            if isinstance(tr, yamlTrace):
                tr = du.yamlTrace2pyrockoTrace(tr)

            lines_dict[source][target] = pltlines.Line2D(tr.get_xdata(),
                    tr.get_ydata())

        return lines_dict 

    def ydata_of_target(self, sources, target):
        if sources == []:
            sources = self.sources
        for source in sources:
            ssmgrm = self.candidates[source][target]
            yield source, ssmgrm.get_xdata(), ssmgrm.get_ydata()

    @property
    def store(self):
        return self.engine.get_store(self.store_id)

    @staticmethod
    def iter_dict(traces_dict, only_values=False):
        """
        Iterate over a 2D-dict, yield each value.
        """
        for key1, key_val in traces_dict.iteritems():
            for key2, val in key_val.iteritems():
                if only_values:
                    yield val
                else:
                    yield key1, key2, val


if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [derec_home + '/fomostos']

    engine = LocalEngine(store_superdirs=store_dirs)
    store_id = 'castor'
    # load stations from file:
    stations = model.load_stations(pjoin(selfdir,
                            '../reference_stations_castor_selection.txt'))

    markers = gui_util.Marker.load_markers(pjoin(selfdir,
                                                '../reference_marker_castor.txt'))

    phase_ids_start = '|'.join(du.get_tabulated_phases(engine,
                                                       store_id, 
                                                       ['p','P']))
    
    # load stations from file:
    # Event==================================================
    event = filter(lambda x: isinstance(x, gui_util.EventMarker), markers)
    assert len(event) == 1
    event = event[0].get_event()
    event.magnitude = 4.3
    event.moment_tensor = moment_tensor.MomentTensor(
                                    m=num.array([[0.0, 0.0, 1.0],
                                                 [0.0, 0.0, 0.0],
                                                 [0.0, 0.0, 0.0]]))


    # generate stations from olat, olon:
    if not stations:
        print 'Generating station distribution.'
        stations = du.station_distribution((event.lat,event.lon),
                                       [[10000., 4], [130000., 8]], 
                                       rotate={3000.:45, 130000.:0})

    targets = du.stations2targets(stations, store_id)

    model = du.get_earthmodel_from_engine(engine, store_id) 

    #TESTSOURCES===============================================
    
    zoffset= 0.
    ref_source = du.event2source(event, 'DC', strike=37.3, dip=30, rake=-3)

    depths=[1800, 2000, 2200]
    #depths=num.linspace(ref_source.depth-zoffset, ref_source.depth+zoffset, 1)

    print depths, '<- depths'

    # Das kann mit als Funktion in TestCaseSetup...
    location_test_sources = [DCSource(lat=ref_source.lat,
                           lon=ref_source.lon,
                           depth=depth,
                           time=event.time,
                           strike=ref_source.strike,
                           dip=ref_source.dip,
                           rake=ref_source.rake,
                           magnitude=event.magnitude) for depth in depths]

    map(lambda x: x.regularize(), location_test_sources)

    reference_request = make_reference_trace(ref_source,
                                                 targets, 
                                                 engine)

    reference_seismograms = du.response_to_dict(reference_request)

    # setup the misfit setup:
    norm = 2.
    taper = trace.CosFader(xfrac=0.2) 
    
    z, p, k = butter(4, (1.*num.pi*2. ,0.4*num.pi*2.) , 
                                       'bandpass', 
                                       analog=True, 
                                       output='zpk')

    #z = num.array(z, dtype=complex)
    z = [complex(zi) for zi in z]
    p = [complex(pi) for pi in p]
    #p = num.array(p, dtype=complex)
    k = complex(k)
    fresponse = trace.PoleZeroResponse(z,p,k)
    fresponse.regularize()

    misfit_setup = trace.MisfitSetup(norm=norm,
                                     taper=taper,
                                     domain='time_domain',
                                     filter=fresponse)

    rise_time=1.
    stf = [[0,rise_time],[0,1]]

    test_case_setup = TestCaseSetup(reference_source=ref_source,
                                    sources=location_test_sources,
                                    targets=targets,
                                    engine=engine, 
                                    store_id=store_id,
                                    misfit_setup=misfit_setup,
                                    source_time_function=stf,
                                    number_of_time_shifts=9,
                                    percentage_of_shift=10.,
                                    phase_ids_start=phase_ids_start) 
    fn = 'sample_test_case_setup.yaml'
    f=open(fn, 'w')

    test_case_setup.regularize()
    test_case_setup.validate()
    f.write(test_case_setup.dump())
    f.close()

    test_case = TestCase( test_case_setup )

    for tr in TestCase.iter_dict(reference_seismograms, only_values=True):
        du.add_random_noise_to_trace(tr, A=0.00001)

    test_case.set_raw_references(reference_seismograms)

    test_case.raw_references = du.apply_stf(test_case.raw_references, 
                            test_case_setup.source_time_function)

    extended_ref_marker = du.chop_ranges(ref_source, 
                                        targets, 
                                        test_case.store,
                                        phase_ids_start,
                                        perc=1.0,
                                        t_shift_frac=0.3)

    test_case.set_reference_markers(extended_ref_marker)

    D = Doer(test_case)

    #test_case.compare_plot( traces_dicts=[test_case.processed_references,
    #                    test_case.processed_candidates],
    #                    focus_first=False)

    #test_case.plot_marker_vs_distance()

    #test_case.plot_z_components(test_case.raw_candidates,
    #        markers_dict=test_case.candidates_markers)
    #
    #test_case.plot_z_components(test_case.raw_references,
    #        sources = [ref_source],
    #        markers_dict=test_case.ref_markers)

    #yaml_trace = yamlTrace()
    #for s,t,tr in TestCase.iter_dict(test_case.raw_references):
    #    yaml_trace.ydata = tr.ydata
    #    yaml_trace.dt = tr.deltat
    #f = open('dump_test.yaml', 'w')
    #f.write(yaml_trace.dump())
    #f.close()
    #plt.show()
    #print 'dumping...'
    test_case.yaml_dump(fn='test_case_dump.yaml')
