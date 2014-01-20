from pyrocko.gf import *
from pyrocko import cake, model, gui_util, util, io, pile, trace, moment_tensor
#from optics import *
from collections import defaultdict

import os
import progressbar
import os
import tempfile
import derec_utils as du
import numpy as num
import copy
import inspect

pjoin = os.path.join

def get_earthmodel_from_engine(engine, store_id):
    return engine.get_store(store_id).config.earthmodel_1d


def stations2targets(stations, store_id):
    '''
    Convert pyrockos original stations into seismosizer targets.
    '''
    targets = []
    for s in stations:
        channels = s.get_channels()
        if channels == []:
            channels = 'NEZ'
        target = [Target(codes=(s.network,s.station,s.location,component),
                                 lat=s.lat,
                                 lon=s.lon,
                                 store_id=store_id,
                                 )for component in channels]
        targets.extend(target)
    return targets


def event2source(event, source_type='MT'):
    '''
    Convert pyrockos original event into seismosizer MT source.

    MT Source magnitude not scaled?!
    returns list of sources
    '''
    if source_type=='MT':
        m = event.moment_tensor._m
        source_event = MTSource(lat=event.lat,
                                   lon=event.lon,
                                   depth=event.depth,
                                   time=event.time,
                                   mnn=float(m[0,0]),
                                   mee=float(m[1,1]),
                                   mdd=float(m[2,2]),
                                   mne=float(m[0,1]),
                                   mnd=float(m[0,2]),
                                   med=float(m[1,2]))
    if source_type=='DC':
        # only one of both possible s,d,r is needed.
        s,d,r = event.moment_tensor.both_strike_dip_rake()[0]
        m = event.moment_tensor.moment_magnitude
        source_event = DCSource(lat=event.lat,
                                lon=event.lon,
                                depth=event.depth,
                                time=event.time,
                                strike=s,
                                dip=d,
                                rake=r,
                                magnitude=m())
    else:
        raise Exception('invalid source type: %s'%source_type)

    return [source_event]


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
    response = engine.process(
            sources=source,
            targets=targets)
    return response.pyrocko_traces()


class Core:
    def __init__(self, markers, stations):
        # Targets================================================
        store_id = 'local1'
        targets = stations2targets(stations, store_id)
        
        # Event==================================================
        event = filter(lambda x: isinstance(x, gui_util.EventMarker), markers)
        assert len(event) == 1
        event = event[0].get_event()
        event.moment_tensor = moment_tensor.MomentTensor(
                                        m=num.array([[1.0, 0.5, 0.0],
                                                     [0.0, 0.5, 0.1],
                                                     [0.0, 0.0, 0.4]]))
        source = event2source(event, 'DC')
        derec_home = os.environ["DEREC_HOME"]
        store_dirs = [derec_home + '/fomostos']

        engine = LocalEngine(store_superdirs=store_dirs)

        reference_seismograms = make_reference_trace(source, targets, engine)
        #TESTSOURCES===============================================
        offset = 0.01
        zoffset= 1000
        print 'z: ',event.depth
        lats=num.arange(event.lat-offset, event.lat+offset, offset/2) 
        print 'lats :',lats 
        lons=num.arange(event.lon-offset, event.lon+offset, offset/1)
        depths=num.arange(event.depth-zoffset, event.depth+zoffset, zoffset/2)
        # only one of both possible s,d,r is needed.
        strike,dip,rake = event.moment_tensor.both_strike_dip_rake()[0]
        m = event.moment_tensor.moment_magnitude
        location_test_sources = [DCSource(lat=lat,
                           lon=lon,
                           depth=depth,
                           time=event.time,
                           strike=strike,
                           dip=dip,
                           rake=rake,
                           magnitude=m()
                           ) for depth in depths for lat in lats for lon in lons]
        print len(location_test_sources)
        #==========================================================
        # Extend P phase markers to 210 p reflection
        #tmpdir = tempfile.mkdtemp(prefix='derec_tmp_', suffix='test')
        primary_phase = cake.PhaseDef('p')

        model = get_earthmodel_from_engine(engine, store_id) 

        #CHOP the reference seimograms once, assuming that markers of P Phases are given:
        extended_ref_marker = du.chop_ranges(model, 
                                             targets,
                                             primary_phase, 
                                             source, 
                                             phase_end=cake.PhaseDef('s'))
        chopped_ref_traces = du.chop_using_markers(reference_seismograms, extended_ref_marker)

        test_case = TestCase(location_test_sources, chopped_ref_traces, targets, engine,mod_parameters=['lat','lon','depth'])
        test_case.request_data()

        # markers hier ueberschreiben. Eigentlich sollen hier die gepicketn Marker verwendet werden. 

        # Das besser in test setup aufrufen:
        #TODO------------------------------------------------------
        # parallelisieren!!!!
        extended_test_marker = du.chop_ranges(model, 
                                              targets,
                                              primary_phase, 
                                              location_test_sources,
                                              phase_end=cake.PhaseDef('P'),
                                              static_offset=8)

        # chop..........................................
        chopped_test_traces = du.chop_using_markers(test_case.response.iter_results(), extended_test_marker) 
 
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
        
        total_misfit = self.calculate_group_misfit(test_case)
        
        test_case.set_misfit(total_misfit)
        print total_misfit

        #memfile = pile.MemTracesFile(parent=None, traces=chopped_test_traces.values()[0])

    def test_vtk(self, test_cases):
        op = OpticBase(test_cases)
        op.numpyrize()

    def calculate_group_misfit(self, test_case):
        # zyklische abhaengigkeit beseitigen!
        candidates = test_case.seismograms
        references = test_case.references
        mfsetups = test_case.misfit_setup
        total_misfit = defaultdict(dict)

        for source in test_case.sources:
            ms = []
            ns = []
            for target in test_case.targets:
                rt = references.values()[0][target]
                
                mf = rt.misfit(candidates=candidates[source].values(), setups=mfsetups)
                for m,n in mf:
                    ms.append(m)
                    ns.append(n)
                print mf
            
            ms = num.array(ms)
            ns = num.array(ns)

            # TODO EXPONENT gleich NORM !!!!!!!!!!
            norm = mfsetups.norm
            M = num.sum(ms**norm)**1/norm
            N = num.sum(ns**norm)**1/norm
                
            total_misfit[source] = M/N

        return total_misfit


class TestCase():
    '''
    In one test case, up to 3 parameters can be modified
    '''
    def __init__(self, sources, references, targets, engine, mod_parameters):
        self.sources = sources
        self.engine = engine
        self.targets = targets
        self.mod_parameters = mod_parameters
        self.references = references

        self.seismograms = {}
        self.results = None
        self.misfit_setup = None
            
    def set_stations(self, stations=[]):
        self.stations = stations 

    def set_misfit_setup(self, setup):
        self.misfit_setup = setup

    def request_data(self):
        print 'requesting data....'
        self.response = self.engine.process(status_callback=self.update_progressbar, sources=self.sources,
                targets=self.targets)
        print 'finished'

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

    def update_progressbar(self, a, b):
        try:
            self.progressbar.update(b)
        except AttributeError:
            self.progressbar = progressbar.ProgressBar(maxval=a).start()
            self.progressbar.update(b)

    def dump_pile(self, fn='test_dumped_seismograms.mseed'):
        pile.make_pile(seismograms.values(), fn=fn)
        

if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    stations = model.load_stations(pjoin(selfdir, '../reference_stations_local.txt'))
    markers = gui_util.Marker.load_markers(pjoin(selfdir, '../reference_marker_local.txt'))
    C = Core(markers=markers, stations=stations)
