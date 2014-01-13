from pyrocko.gf.seismosizer import *
from pyrocko import cake, model, gui_util, util, io, pile, trace, moment_tensor
import os
import tempfile
import derec_utils as du
import numpy as num
import copy
import inspect

pjoin = os.path.join

def equal_attributes(o1, o2):
    return o1.__dict__ == o2.__dict__

def set_refine_parameter(ref_event, **kwargs):
    events = {}
    for k, vals in kwargs.iteritems():
        for val in vals:
            event_copy = copy.copy(ref_event)
            exec('event_copy.%s=%s' % (k, val))
            events[val]=event_copy

    return events


def generate_test_events(event, **kwargs):
    events = {}
    for key in kwargs:
        if not key in inspect.getargspec(model.Event.__init__).args:
            raise Exception(''''key %s not possible to refine. Possible options 
                    are: %s''' % (key, inspect.getargspec(model.Event.__init__).args))
        if len(kwargs.items())>=2:
            raise Exception('''Too many refine parameters. Give one parameter,
                                one range''')
        else:
            events = set_refine_parameter(event, **kwargs)       
    return events


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
        pass

        store_id = 'crust2_dd'

        event = filter(lambda x: isinstance(x, gui_util.EventMarker), markers)
        assert len(event) == 1
        event = event[0].get_event()

        event.moment_tensor = moment_tensor.MomentTensor(m=num.array([[1.0, 0.5, 0.0],
                                                                      [0.0, 0.5, 0.1],
                                                                      [0.0, 0.0, 0.4]]))
        
        reference_seismograms = make_reference_trace(event,
                                                     stations,
                                                     store_id='crust2_dd')

        tmpdir = tempfile.mkdtemp(prefix='derec_tmp_', suffix='test')    
        io.save(reference_seismograms, filename_template=pjoin(tmpdir, 'ref.mseed'))

        test_depths = num.arange(1000,8000,2000)
        test_events = generate_test_events(event, depth=test_depths)
        test_case = TestCase(test_events, stations)
        test_case.request_data()
        test_seismograms = test_case.get_seismograms()

        # Request earthmodel
        model_request = request_earthmodel(test_case.store_id)
        model = model_request.store_configs[0].earthmodel_1d

        map(lambda s: s.set_event_relative_data(event), stations)

        # Extend P phase markers to 210 p reflection
        #latest_phase = cake.PhaseDef('pPv210p')
        #latest_phase = cake.PhaseDef('S')
        primary_phase = cake.PhaseDef('P')
        #primary_phase = cake.PhaseDef('p')

        # markers hier ueberschreiben. Eigentlich sollen hier die gepicketn Marker verwendet werden. 

        # Das besser in test setup aufrufen:
        extended_test_marker = du.chop_ranges(model, 
                                              stations,
                                              primary_phase, 
                                              test_events)

        # chop..........................................
        equal_test_event = filter( lambda e: equal_attributes(e, event), test_events.values())[0]
        chopped_ref_traces = []
        chopped_ref_traces.extend(du.chop_using_markers(reference_seismograms, extended_test_marker[equal_test_event]))
        
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
        
        test_case.set_results(total_misfit)
       
        du.plot_misfit_dict(total_misfit)

        #memfile = pile.MemTracesFile(parent=None, traces=chopped_test_traces.values()[0])

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
        self.keys = {}
        self.seismograms = {}
        self.results = None
        self.misfit_setup = None
            
    def set_stations(self, stations=[]):
        self.stations = stations 

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

    def set_seismograms(self, seismograms):
        self.seismograms = seismograms

    def set_results(self, results):
        self.results = results

    def dump_requests(self):
        for r in self.requests:
            fn = 'test.yaml'
            f = open(fn, 'w')
            f.write(r.dump())
            f.close()

    def dump_pile(self, fn='test_dumped_seismograms.mseed'):
        pile.make_pile(seismograms.values(), fn=fn)
        


if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    stations = model.load_stations(pjoin(selfdir, '../reference_stations.txt'))
    markers = gui_util.Marker.load_markers(pjoin(selfdir, '../reference_marker.txt'))
    C = Core(markers=markers, stations=stations)
