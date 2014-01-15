from pyrocko.gf import *
from pyrocko import cake, model, gui_util, util, io, pile, trace, moment_tensor
#from optics import *

import os
import tempfile
import derec_utils as du
import numpy as num
import copy
import inspect

pjoin = os.path.join


def stations2targets(stations, store_id):
    '''
    Convert pyrockos original stations into seismosizer targets.
    '''
    targets = []
    for s in stations:
        targets.append(gf.Target(codes=(s.network,s.station,s.location,component),
                                 lat=s.latitude,
                                 lon=s.longitude,
                                 store_id=store_id,
                                 )for component in s.channels)
    return targets


def event2source(event, store_id):
    '''
    Convert pyrockos original event into seismosizer source.
    '''
    m = event.moment_tensor._m
    source_event = gf.MTSource(lat=event.latitude,
                               lon=event.longitude,
                               depth=event.depth,
                               magnitude=event.magnitude,
                               mnn=m[0,0],
                               mee=m[1,1],
                               mdd=m[2,2],
                               mne=m[0,1],
                               mnd=m[0,2],
                               med=m[1,2])
    return source_event


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


def generate_test_events(event, **kwargs):
    '''
    Returns dict with keys=z , values= event objects
    '''
    for key,arg in kwargs.iteritems():
        
    #events = {}
    #for key in kwargs:
    #    if not key in inspect.getargspec(model.Event.__init__).args:
    #        raise Exception(''''key %s not possible to refine. Possible options 
    #                are: %s''' % (key, inspect.getargspec(model.Event.__init__).args))
    #    if len(kwargs.items())>=2:
    #        raise Exception('''Too many refine parameters. Give one parameter,
    #                            one range''')
    #    else:
    #        events = set_refine_parameter(event, **kwargs)       
    #return events


def make_reference_trace(event, stations, store_id='crust2_dd'):
    reference_seismograms = []
    m = event.moment_tensor._m
    for stat in stations:
        reference_seis_req = SeismosizerRequest(store_id=store_id,
                                                source_lat=event.lat,
                                                source_lon=event.lon,
                                                source_depth=event.depth,
                                                receiver_lat=stat.lat,
                                                receiver_lon=stat.lon,
                                                source_time=event.time,
                                                net_code=stat.network,
                                                sta_code=stat.station,
                                                loc_code=stat.location,
                                                mnn=m[0,0],
                                                mee=m[1,1],
                                                mdd=m[2,2],
                                                mne=m[0,1],
                                                mnd=m[0,2],
                                                med=m[1,2])

        reference_seismograms.extend(request_seismogram(reference_seis_req).traces)
    return reference_seismograms


class Processor():
    def __init__(self, test_case):
        pass


class Core:
    def __init__(self, markers, stations):
        # Targets================================================
        store_id = 'crust2_dd'
        targets = stations2targets(stations, store_id)
        
        # Event==================================================
        event = filter(lambda x: isinstance(x, gui_util.EventMarker), markers)
        assert len(event) == 1
        event = event[0].get_event()
        event.moment_tensor = moment_tensor.MomentTensor(m=num.array([[1.0, 0.5, 0.0],
                                                                      [0.0, 0.5, 0.1],
                                                                      [0.0, 0.0, 0.4]]))

        reference_seismograms = make_reference_trace(source,
                                                     stations,
                                                     store_id='crust2_dd')

        tmpdir = tempfile.mkdtemp(prefix='derec_tmp_', suffix='test')    

        # Extend P phase markers to 210 p reflection
        primary_phase = cake.PhaseDef('P')

        # Request earthmodel
        model_request = request_earthmodel(test_case.store_id)
        model = model_request.store_configs[0].earthmodel_1d

        #CHOP the reference seimograms once, assuming that markers of P Phases are given:
        map(lambda s: s.set_event_relative_data(event), stations)
        extended_ref_marker = du.chop_ranges(model, 
                                              stations,
                                              primary_phase, 
                                              event)

        chopped_ref_traces = []
        chopped_ref_traces.extend(du.chop_using_markers(reference_seismograms, extended_ref_marker[event]))

        io.save(reference_seismograms, filename_template=pjoin(tmpdir, 'ref.mseed'))

        test_depths = num.arange(1000,8000,2000)
        test_lats = num.arange(event.latitude-1, event.latitude+1, 0.4)

        event_copy = event.copy()            
        event_copy.latitude = lat
        test_events = generate_test_events(event_copy, depth=test_depths)
        test_case = TestCase(test_events, stations)
        test_case.request_data()
        test_seismograms = test_case.get_seismograms()

        map(lambda s: s.set_event_relative_data(event), stations)

        # markers hier ueberschreiben. Eigentlich sollen hier die gepicketn Marker verwendet werden. 

        # Das besser in test setup aufrufen:
        extended_test_marker = du.chop_ranges(model, 
                                              stations,
                                              primary_phase, 
                                              test_events)

        # chop..........................................
        
        chopped_test_traces = {}
        for e in test_events.values():
            chopped_test_traces[e] = du.chop_using_markers(test_seismograms[e], extended_test_marker[e]) 
        test_case.set_seismograms(chopped_test_traces)
 
        # Misfit.........................................
        norm = 2
        taper = trace.CosFader(xfade=3)
        fresponse = trace.FrequencyResponse()
        setup = trace.MisfitSetup(norm=norm,
                                  taper=taper,
                                  domain='time_domain',
                                  freqlimits=(1,2,20,40),
                                  frequency_response=fresponse)
        
        test_case.set_misfit_setup(setup)
        
        total_misfit = self.calculate_group_misfit(chopped_ref_traces,
                                                   test_case)
        
        test_case.set_misfit(total_misfit)

        #memfile = pile.MemTracesFile(parent=None, traces=chopped_test_traces.values()[0])

    def test_vtk(self, test_cases):
        op = OpticBase(test_cases)
        op.numpyrize()

    def calculate_group_misfit(self, traces, test_case):
        candidates = test_case.seismograms
        mfsetups = test_case.misfit_setup
        total_misfit = {}
        for d,tts in candidates.items():
            ms = []
            ns = []
            for rt in traces:
                for tt in tts:
                    if rt.nslc_id==tt.nslc_id:
                        # TODO: nach candidates vorsortieren.
                        mf = rt.misfit(candidates=[tt], setups=mfsetups)
                        for m,n in mf:
                            ms.append(m)
                            ns.append(n)
            
            ms = num.array(ms)
            ns = num.array(ns)

            # TODO EXPONENT gleich NORM !!!!!!!!!!
            M = num.sqrt(num.sum(ms**2))
            N = num.sqrt(num.sum(ns**2))
                
            total_misfit[d] = M/N

        return total_misfit


class TestCase():
    def __init__(self, events, stations, store_id='crust2_dd'):
        
        self.stations = stations 
        self.events = events
        self.store_id = store_id
        self.key = ''
        self.seismograms = {}
        self.results = None
        self.misfit_setup = None
            
    def set_stations(self, stations=[]):
        self.stations = stations 

    def set_key(self):
        '''key is the parameter, that is varied'''
        self.key = key

    def set_misfit_setup(self, setup):
        self.misfit_setup = setup

    def request_data(self):
        for k, event in self.events.iteritems():
            print 'requesting data for', event
            seismos= []
            for stat in self.stations:
                s = du.request_data(stat, event, self.store_id)
                seismos.extend(s)

            self.seismograms[event] = seismos

    def get_seismograms(self):
        return self.seismograms

    def get_misfit_array(self):
        '''
        Should return a numpy array containing the results, after these have 
        been sorted by the varying parameter (key).
        '''
        misfit_array = num.zeros(len(self.results))
        return num.array([sorted(self.results.keys(), key=operator.attrgetter(self.key))])

    def set_seismograms(self, seismograms):
        self.seismograms = seismograms

    def set_misfit(self, results):
        self.results = results

    def dump_requests(self):
        for r in self.requests:
            fn = 'test.yaml'
            f = open(fn, 'w')
            f.write(r.dump())
            f.close()

    def dump_pile(self, fn='test_dumped_seismograms.mseed'):
        pile.make_pile(seismograms.values(), fn=fn)
        
#class TestCase3D(TestCase):
#    def __init__(self, events, stations, store_id='')
#        self.super(TestCase).__init__(events, stations, store_id)
    

if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    stations = model.load_stations(pjoin(selfdir, '../reference_stations.txt'))
    markers = gui_util.Marker.load_markers(pjoin(selfdir, '../reference_marker.txt'))
    C = Core(markers=markers, stations=stations)
