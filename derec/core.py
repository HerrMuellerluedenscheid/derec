from pyrocko.gf import *
from pyrocko import cake, model, gui_util, util, io, pile, trace, moment_tensor
from pyrocko import orthodrome
from vtkOptics import *
from collections import defaultdict
from matplotlib import cm
from gmtpy import griddata_auto

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

    elif source_type=='DC':
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
                                magnitude=5)

    elif source_type=='EX':
        m = event.moment_tensor.moment_magnitude
        source_event = ExplosionSource(lat=event.lat,
                                       lon=event.lon,
                                       depth=event.depth,
                                       time=event.time,
                                       magnitude=5)
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
    

def make_reference_markers(source, targets, model):
    
    assert len(source) == 1

    ref_marker = defaultdict(dict)
    phases_start = ['p','P']
    phases_start = [cake.PhaseDef(ph) for ph in phases_start]

    phases_end = ['s', 'S']
    phases_end = [cake.PhaseDef(ph) for ph in phases_end]
    
    for s in source:
        for target in targets:
            dist = orthodrome.distance_accurate50m(s, target)*cake.m2d
            tmin = min(model.arrivals([dist], 
                                    phases_start,
                                    zstart=s.depth,
                                    zstop=s.depth), key=lambda x: x.t).t

            tmax = min(model.arrivals([dist], 
                                    phases_end, 
                                    zstart=s.depth,
                                    zstop=s.depth), key=lambda x: x.t).t

            tmin += s.time
            tmax += s.time
            assert tmin!=tmax
            
            m = gui_util.PhaseMarker(nslc_ids=target.codes, 
                                    tmin=tmin,
                                    tmax=tmax,
                                    kind=1,
                                    event=source,
                                    phasename='range')

            ref_marker[s][target] = m
    return ref_marker


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
                                        m=num.array([[1.0, 0.0, 0.0],
                                                     [0.0, 1.0, 0.0],
                                                     [0.0, 0.0, 1.]]))
    
        source = list(event2source(event, 'DC'))

        derec_home = os.environ["DEREC_HOME"]
        store_dirs = [derec_home + '/fomostos']

        engine = LocalEngine(store_superdirs=store_dirs)
        model = get_earthmodel_from_engine(engine, store_id) 

        extended_ref_marker = make_reference_markers(source, targets, model)
        reference_seismograms = make_reference_trace(source, targets, engine)
        
        #TESTSOURCES===============================================
        offset = 0.1
        zoffset= 10000.

        lats=num.arange(event.lat-0.5*offset, event.lat+0.5*offset, offset/1) 
        lons=num.arange(event.lon-0.5*offset, event.lon+0.5*offset, offset/1)
        #lons=[event.lon]
        #lats=[event.lat]
        #depths = [event.depth]
        
        print lats, lons, '<- lats, lons'

        depths=num.arange(event.depth-zoffset, event.depth+zoffset, zoffset/3)
        strike,dip,rake = event.moment_tensor.both_strike_dip_rake()[1]
        m = event.moment_tensor.moment_magnitude
        location_test_sources = [DCSource(lat=lat,
                               lon=lon,
                               depth=depth,
                               time=event.time,
                               strike=strike,
                               dip=dip,
                               rake=rake,
                               magnitude=5) for depth in depths for lat in lats for lon in lons]
        
        #==========================================================

        test_case = TestCase(location_test_sources, 
                             targets, 
                             engine, 
                             store_id, 
                             test_parameters={'depth':depths, 'lat':lats, 'lon':lons})

        test_case.request_data()
        extended_test_marker = du.chop_ranges(test_case, 'p|P|Pv_35p', 's|S|Sv_35s', t_start_shift=-2, t_end_shift=-2)
        test_case.references = du.chop_using_markers(reference_seismograms, extended_ref_marker)

        # chop..........................................
        test_case.seismograms = du.chop_using_markers(test_case.response.iter_results(), extended_test_marker) 

        # Misfit.........................................
        norm = 2.
        taper = trace.CosFader(xfade=3) # Seconds or samples?
        fresponse = trace.FrequencyResponse()
        setup = trace.MisfitSetup(norm=norm,
                                  taper=taper,
                                  domain='time_domain',
                                  freqlimits=(2,4,20,40),
                                  filter=fresponse)
        
        test_case.set_misfit_setup(setup)
        total_misfit = self.calculate_group_misfit(test_case)
        test_case.set_misfit(total_misfit)
        #test_case.plot1d()
        test_case.contourf()

        test_tin = TestTin([test_case])
        optics = OpticBase(test_tin)
        #optics.plot_1d(fix_parameters={'lat':event.lat, 'lon':event.lon})

    def calculate_group_misfit(self, test_case):
        # zyklische abhaengigkeit beseitigen!
        candidates = test_case.seismograms
        references = test_case.references
        assert len(references.items())==1

        mfsetups = test_case.misfit_setup
        total_misfit = defaultdict(dict)

        for source in test_case.sources:
            ms = []
            ns = []
            
            #print 'new_source'
            for target in test_case.targets:
                rt = references.values()[0][target]
                mf = rt.misfit(candidates=[candidates[source][target]], setups=mfsetups)
                for m,n in mf:
                    ms.append(m)
                    ns.append(n)
                
            ms = num.array(ms)
            ns = num.array(ns)

            # TODO EXPONENT gleich NORM !!!!!!!!!!
            norm = mfsetups.norm
            M = num.sum(ms)
            N = num.sum(ns)
            #M = num.power(num.sum(num.power(ms, norm)), 1./norm)
            #N = num.power(num.sum(num.power(ns, norm)), 1./norm)
                
            total_misfit[source] = M/N

        return total_misfit


class TestTin():
    def __init__(self, test_cases=[]):
        self.assertSameParameters(test_cases)
        self.test_cases = test_cases
        self.test_parameters = self.test_cases[0].test_parameters
        self.update_dimensions(self.test_cases)
        self.x_range = None
        self.y_range = None
        self.z_range = None

    def add_test_case(test_case):
        self.assertSameParameters(test_case)
        self.test_cases.extend(test_case)
        self.update_dimensions([test_case])

    def numpyrize_1d(self, fix_parameters={}):
        '''Make 1dimensional numpy array

        fix_parameters is a dict {parameter1:value, parameter2:value}

        ..examples:
            numpyrize_1d({latitude:10, depth:1000})

        '''
        assert all(k in self.test_parameters for k in fix_parameters.keys())
        if not len(fix_parameters.keys())==len(self.test_parameters)-1:
            raise Exception('Expected %s fix_parameters, got %s' % (len(self.test_parameters)-1, len(fix_parameters.keys())))

        x = []
        y = []
        x_key = ''.join(set(self.test_parameters) - set(fix_parameters.keys()))
        
        for tc in self.test_cases:
            misfits = tc.misfits
            for source in tc.sources:
                if all(getattr(source, k)==v for k,v in fix_parameters.items()):
                    x.append(getattr(source, x_key))
                    y.append(misfits[source])
                else:
                    continue
        return num.array(x), num.array(y)

    def numpyrize_2d(self, fix_parameters={}):
        '''Make 1dimensional numpy array

        fix_parameters is a dict {parameter1:value, parameter2:value}

        ..examples:
            numpyrize_2d({latitude:10, depth:1000})

        '''
        assert all(k in self.test_parameters for k in fix_parameters.keys())
        if not len(fix_parameters.keys())==len(self.test_parameters)-1:
            raise Exception('Expected %s fix_parameters, got %s' % (len(self.test_parameters)-1, len(fix_parameters.keys())))

        data = num.ndarray(shape=(len(tc),len()))

        x = []
        y = []
        z = []

        # TODO: Ordering of keys needs revision.
        x_key, y_key = set(self.test_parameters) - set(fix_parameters.keys())
        
        for tc in self.test_cases:
            misfits = tc.misfits
            for source in tc.sources:
                if all(getattr(source, k)==v for k,v in fix_parameters.items()):
                    D[getattr(source, xkey)][getattr(source,ykey)]=misfits[source]

        return x, y

    def update_dimensions(self, test_cases):
        for tc in test_cases:
            self.x_range = num.union1d(self.x_range, tc.test_parameters.items()[0].values())
            self.y_range = num.union1d(self.y_range, tc.test_parameters.items()[1].values())
            self.z_range = num.union1d(self.z_range, tc.test_parameters.items()[2].values())


    def assertSameParameters(self, test_cases):
        if not isinstance(test_cases, list):
            test_cases = list(test_cases)

        assert all(set(x.test_parameters)==set(test_cases[0].test_parameters) for x in test_cases)



class TestCase():
    '''
    In one test case, up to 3 parameters can be modified
    '''
    def __init__(self, sources, targets, engine, store_id, test_parameters):
        self.test_grid = None
        self.sources = sources
        self.engine = engine
        self.targets = targets
        self.test_parameters = test_parameters 
        self.store_id = store_id

        self.seismograms = {}
        self.misfits = None
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

    def set_seismograms(self, seismograms):
        self.seismograms = seismograms

    def set_misfit(self, misfits):
        self.misfits=misfits 

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

    def snuffle(self):
        trace.snuffle(self.seismograms)

    def numpy_it(self, **kwargs):
        '''
        Irgendwann muss ein ordering mit angegeben werden. Dicts sind ungeordnet.
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

    def plot1d(self):
        self.numpy_it()
        plt.plot(self.num_array[1], self.num_array[3])
        plt.show()

    def contourf(self):
        
        import pdb
        pdb.set_trace()
        self.numpy_it(order=['lat', 'lon','depth'])

        x=self.num_array[0]
        y=self.num_array[1]
        z=self.num_array[2]
        v=self.num_array[3]

        X,Y,Z = griddata_auto(x,y,v)
        plt.contourf(X,Y, Z, cmap=cm.bone_r)
        plt.xlabel(self.xkey)
        plt.ylabel(self.ykey)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('L%s misfit'%int(self.misfit_setup.norm))

        plt.show()

    @property
    def store(self):
        return self.engine.get_store(self.store_id)


if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    stations = model.load_stations(pjoin(selfdir, '../reference_stations_local.txt'))
    markers = gui_util.Marker.load_markers(pjoin(selfdir, '../reference_marker_local.txt'))
    C = Core(markers=markers, stations=stations)
