from pyrocko.gf.seismosizer import *
from pyrocko import cake, model, gui_util, util, io, pile, trace, moment_tensor
import os
import tempfile
import derec_utils as du
import numpy as num
import copy

pjoin = os.path.join

def generate_test_events(event, **kwargs):
    events = {}
    assert len(kwargs.items())==1
    for k, vals in kwargs.items():
        assert isinstance(vals, num.ndarray)
        for val in vals:
            event_copy = copy.copy(event)
            setattr(event_copy, k, val)
            events[val] = event_copy
    
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


class Core:
    def __init__(self, markers, stations):

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
        # collect all test cases
        test_depths = num.arange(1000,8000,1000)
        test_events = generate_test_events(event, depth=test_depths)
        test_case = TestCaseBase(test_events, stations)

        # erzeuge test events:
        #test_case.set_refine_parameter(source_depth=test_depths)
        test_case.request_data()

        test_seismograms = test_case.get_seismograms()

        # Check if ydata in correct depths are equal
        #for ts in test_seismograms[event.depth]:
        #    for rs in reference_seismograms:
        #        if ts.nslc_id==rs.nslc_id:
        #            assert (ts.get_ydata()==rs.get_ydata()).all()

        # Request earthmodel
        model_request = request_earthmodel(store_id)
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
                                              event, 
                                              primary_phase, 
                                              test_depths)

        #extended_markers = list(du.extend_phase_markers(markers=markers,
        #                                                phase=latest_phase,
        #                                                stations=stations,
        #                                                event=event, model=model))

        # chop
        chopped_ref_traces = []
        chopped_ref_traces.extend(du.chop_using_markers(reference_seismograms, extended_test_marker[event.depth]))
        
        chopped_test_traces = {}
        for d in test_depths:
            chopped_test_traces[d] = du.chop_using_markers(test_seismograms[d], extended_test_marker[d]) 

        norm = 2
        taper = trace.CosFader(xfade=3)
        fresponse = trace.FrequencyResponse()
        setup = trace.MisfitSetup(norm=norm,
                                  taper=taper,
                                  domain='time_domain',
                                  freqlimits=(1,2,20,40),
                                  frequency_response=fresponse)
        
        
        total_misfit = self.calculate_group_misfit(chopped_ref_traces,
                                     chopped_test_traces,
                                     setup) 
       
        du.plot_misfit_dict(total_misfit)

        # testweise nur element 0
        #memfile = pile.MemTracesFile(parent=None, traces=chopped_test_traces.values()[0])
        #p = pile.Pile()
        #inj = pile.Injector(p)
        #seismograms[0].snuffle()
        #inj.inject(seismograms[0])
        #from pyrocko.snuffler import snuffle
        #snuffle(memfile)

    def calculate_group_misfit(self, traces, candidates, mfsetups):
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

            M = num.sqrt(num.sum(ms**2))
            N = num.sqrt(num.sum(ns**2))
                
            total_misfit[d] = M/N

        return total_misfit
            

class TestCaseBase():
    def __init__(self, events, stations, store_id='crust2_dd'):
        
        self.stations = stations 
        self.events = events
        self.store_id = store_id
        #self.refine_parameter_options = dir(self.events[0])
        self.keys = {}
        self.seismograms = {}
            
    def set_stations(self, stations=[]):
        self.stations = stations 

    def set_refine_parameter(self, **kwargs): 
        for key in kwargs:
            if not key in dir(SeismosizerRequest):
                raise Exception('key %s not possible to refine' % key)
            else:
                self.keys[key] = kwargs[key]

    def request_data(self):
        for event in self.events:
            seismos= []
            for stat in self.stations:
                s = du.request_data(stat, event, self.store_id)
                seismos.extend(s)

            self.seismograms[event] = seismos



        '''
        self.traces = None
        event = self.event
        store_id = self.store_id
        source_lat = event.lat
        source_lon = event.lon
        source_depth = event.depth
        source_time = event.time

        mt = event.moment_tensor
        mnn = mt._m[0,0]
        mee = mt._m[1,1]
        mdd = mt._m[2,2] 
        mne = mt._m[0,1]
        mnd = mt._m[0,2]
        med = mt._m[1,2]
        
        for k,values in self.keys.items():
            for value in values:
                exec('%s=%s'%(k, value))
                for stat in self.stations:
                    test_seis_req = SeismosizerRequest(store_id=store_id,
                                                       source_lat=source_lat,
                                                       source_lon=source_lon,
                                                       source_depth=source_depth,
                                                       receiver_lat=stat.lat,
                                                       receiver_lon=stat.lon,
                                                       source_time=source_time,
                                                       net_code=stat.network,
                                                       sta_code=stat.station,
                                                       loc_code=stat.location,
                                                       mnn=mnn,
                                                       mee=mee,
                                                       mdd=mdd,
                                                       mne=mne,
                                                       mnd=mnd,
                                                       med=med)
                    
                    #try:
                    #    self.seismograms[stat.nsl()]
                    #except KeyError:
                    #    self.seismograms[stat.nsl()] = []
                    try:
                        self.seismograms[value]
                    except KeyError:
                        self.seismograms[value] = []

                    self.seismograms[value].extend(request_seismogram(test_seis_req).traces[:])
    '''
    def get_seismograms(self):
        return self.seismograms

    def dump_requests(self):
        for r in self.requests:
            fn = 'test.yaml'
            f = open(fn, 'w')
            f.write(r.dump())
            f.close()


if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    stations = model.load_stations(pjoin(selfdir, '../reference_stations.txt'))
    markers = gui_util.Marker.load_markers(pjoin(selfdir, '../reference_marker.txt'))
    C = Core(markers=markers, stations=stations)
