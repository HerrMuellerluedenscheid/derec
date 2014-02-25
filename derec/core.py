from pyrocko.gf import *
from pyrocko import model, gui_util, pile, trace, moment_tensor
from vtkOptics import *
from collections import defaultdict
from matplotlib import cm
from gmtpy import griddata_auto
from scipy.signal import butter
from guts import *

import time
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
import progressbar
import os
import derec_utils as du
import numpy as num
import copy

pjoin = os.path.join
km = 1000.

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
    
    map(lambda x: x.regularize(), targets)
    return targets


def event2source(event, source_type='MT', rel_north_shift=0., rel_east_shift=0.):
    '''
    Convert pyrockos original event into seismosizer MT source.

    MT Source magnitude not scaled?!
    returns list of sources
    '''
    rel_n_deg, rel_e_deg = du.lat_lon_relative_shift(event.lat, event.lon,
                                rel_north_shift, rel_east_shift)

    if source_type=='MT':
        m = event.moment_tensor._m
        source_event = MTSource(lat=rel_n_deg,
                                   lon=rel_e_deg,
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
        source_event = DCSource(lat=rel_n_deg,
                                lon=rel_e_deg,
                                depth=event.depth,
                                time=event.time,
                                strike=s,
                                dip=d,
                                rake=r,
                                magnitude=event.magnitude)

    elif source_type=='EX':
        m = event.moment_tensor.moment_magnitude
        source_event = ExplosionSource(lat=rel_n_deg,
                                       lon=rel_e_deg,
                                       depth=event.depth,
                                       time=event.time,
                                       magnitude=event.magnitude)
    else:
        raise Exception('invalid source type: %s'%source_type)

    
    source_event.regularize()
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
    return response
    

class Core:
    def __init__(self, markers, stations=None):
        # Targets================================================
        store_id = 'castor'

        if store_id=='local1':
            phase_ids_start = 'p|P|Pv20p|Pv35p'
            phase_ids_end =   's|S|Sv20s|Sv35s'

        if store_id=='very_local':
            phase_ids_start = 'p|P|Pv3p|Pv8p|Pv20p|Pv35p'
            phase_ids_end =   's|S|Sv3s|Sv8s|Sv20s|Sv35s'

        if store_id=='very_local_20Hz':
            phase_ids_start = 'begin_fallback|p|P|Pv1p|Pv3p|Pv8p|Pv20p|Pv35p'
            phase_ids_end =   's|S|Sv1s|Sv3s|Sv8s|Sv20s|Sv35s'

        if store_id=='very_local_20Hz':
            phase_ids_start = 'begin_fallback|p|P|Pv1p|Pv3p|Pv8p|Pv20p|Pv35p'
            phase_ids_end =   's|S|Sv1s|Sv3s|Sv8s|Sv20s|Sv35s'

        if store_id=='castor':
            # bug?! bei Pv1.5p gibt's bei nahen Entfernungen ein index ot of
            # bounds
            phase_ids_start = 'p|P|Pv12.5p|Pv2.5p|Pv18.5p|Pv20p|Pv35p'
            phase_ids_end= 's|S|Sv12.5s|Sv2.5s|Sv18.5s|Sv20s|Sv35s'

        # Event==================================================
        event = filter(lambda x: isinstance(x, gui_util.EventMarker), markers)
        assert len(event) == 1
        event = event[0].get_event()
        event.magnitude = 4.3
        event.moment_tensor = moment_tensor.MomentTensor(
                                        m=num.array([[1.0, 0.0, -0.0],
                                                     [0.0, 1.0, 0.0],
                                                     [0.0, 0.0, 1.0]]))
    

        # generate stations from olat, olon:
        if not stations:
            stations = du.station_distribution((event.lat,event.lon),
                                           [[10000., 4], [130000., 8]], 
                                           rotate={3000.:45, 130000.:0})

        targets = stations2targets(stations, store_id)

        derec_home = os.environ["DEREC_HOME"]
        store_dirs = [derec_home + '/fomostos']

        engine = LocalEngine(store_superdirs=store_dirs)
        model = get_earthmodel_from_engine(engine, store_id) 

        
        #TESTSOURCES===============================================
        
        offset = 0.05

        #offset von km in degrees umrechnen
        zoffset= 1000.

        lats=num.arange(event.lat-offset, event.lat+offset, offset/5) 
        lons=num.arange(event.lon-offset, event.lon+offset, offset/5)
        
        #lons = [event.lon]
        #lats = [event.lat]
        #depths=num.arange(event.depth-zoffset, event.depth+zoffset, zoffset/5)
        depths = [event.depth]
        print lats, '<- lats'
        print lons, '<- lons'
        print depths, '<- depths'
        print event.lat, '<- event lat'
        print event.lon, '<- event lon'
        print event.depth, '<- event depth'

        strike,dip,rake = event.moment_tensor.both_strike_dip_rake()[0]
        print strike, dip, rake
        location_test_sources = [DCSource(lat=lat,
                               lon=lon,
                               depth=depth,
                               time=event.time,
                               strike=strike,
                               dip=dip,
                               rake=rake,
                               magnitude=event.magnitude) for depth in depths 
                                        for lat in lats for lon in lons]

        for s in location_test_sources:
            s.regularize()
        #==========================================================

        test_case = TestCase(location_test_sources, 
                             targets, 
                             engine, 
                             store_id, 
                             test_parameters={'depth':depths, 
                                              'lat':lats, 
                                              'lon':lons})

        test_case.request_data()

        test_case.ref_source = list(event2source(event, 'DC', rel_north_shift=0*km,
                                                    rel_east_shift=0*km))

        print id(test_case.targets[0])
        reference_seismograms = make_reference_trace(test_case.ref_source,
                                                     test_case.targets, 
                                                     engine)
        
        print id(test_case.targets[0])
        extended_ref_marker = du.chop_ranges(test_case.ref_source, 
                                             test_case.targets, 
                                             test_case.store,
                                             phase_ids_start,
                                             phase_ids_end,
                                             t_shift_frac=0.10)


        print('test data marker....')
        t1 = time.time()
        extended_test_marker = du.chop_ranges(test_case.sources,
                                              test_case.targets,
                                              test_case.store,
                                              phase_ids_start, 
                                              phase_ids_end, 
                                              t_shift_frac=0.10)
        
        print time.time()-t1
        test_case.test_markers = extended_test_marker
        test_case.ref_markers = extended_ref_marker

        print('chopping ref....')
        print id(test_case.targets[0])
        test_case.references = du.chop_using_markers(reference_seismograms.iter_results(),
                                                     extended_ref_marker)

        print('chopping cand....')
        test_case.seismograms = du.chop_using_markers(test_case.response.iter_results(),
                                                     extended_test_marker) 

        norm = 2.
        #taper = trace.CosFader(xfade=4) # Seconds or samples?
        taper = trace.CosFader(xfrac=0.1) 
        
        z, p, k = butter(4, (2.*num.pi*2. ,0.4*num.pi*2.) , 
                           'bandpass', 
                           analog=True, 
                           output='zpk')

        z = num.array(z, dtype=complex)
        p = num.array(p, dtype=complex)
        k = num.complex(k)
        fresponse = trace.PoleZeroResponse(z,p,k)
        fresponse.regularize()
        setup = trace.MisfitSetup(norm=norm,
                                  taper=taper,
                                  domain='time_domain',
                                  filter=fresponse)

        test_case.set_misfit_setup(setup)
        du.calculate_misfit(test_case)
        #test_case.yaml_dump()

        # Display results===================================================
        #test_case.plot1d(order, event.lon)
        test_case.contourf({'depth':event.depth})
        #test_case.check_plot({'lat':event.lat, 'depth':event.depth})

        #test_tin = TestTin([test_case])
        #optics = OpticBase(test_tin)
        #optics.plot_1d(fix_parameters={'lat':event.lat, 'lon':event.lon})


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
        assert all([k in self.test_parameters for k in fix_parameters.keys()])
        if not len(fix_parameters.keys())==len(self.test_parameters)-1:
            raise Exception('Expected %s fix_parameters, got %s' % 
                        (len(self.test_parameters)-1, len(fix_parameters.keys())))

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
            raise Exception('Expected %s fix_parameters, got %s' % 
                        (len(self.test_parameters)-1, len(fix_parameters.keys())))

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
            self.x_range = num.union1d(self.x_range, 
                            tc.test_parameters.items()[0].values())
            self.y_range = num.union1d(self.y_range, 
                            tc.test_parameters.items()[1].values())
            self.z_range = num.union1d(self.z_range, 
                            tc.test_parameters.items()[2].values())


    def assertSameParameters(self, test_cases):
        if not isinstance(test_cases, list):
            test_cases = list(test_cases)

        assert all(set(x.test_parameters)==set(test_cases[0].test_parameters) 
                            for x in test_cases)


class TestCase(Object):
    '''
    In one test case, up to 3 parameters can be modified
    '''

    
    sources = List.T(Source.T()) 
    targets = List.T(Target.T()) 
    engine = Engine.T()
    store_id = String.T()
    test_parameters = List.T(String.T())
    misfit_setup = trace.MisfitSetup.T()

    def __init__(self, sources=sources, targets=targets, engine=engine, 
                        store_id=store_id, test_parameters=test_parameters):
        self.targets=targets
        # sollte unnoetig sein:
        self.sources =sources
        self.engine = engine
        self.test_parameters = test_parameters 
        self.store_id = store_id

        self.processed_references = defaultdict(dict)
        self.references = {}
        self.processed_candidates = defaultdict(dict)
        self.seismograms = {}
        self.misfits = None
        self.misfit_setup = None
            
    def set_stations(self, stations=[]):
        self.stations = stations 

    def set_misfit_setup(self, setup):
        self.misfit_setup = setup
        #self.misfit_setup.regularize(depth=10)

    def request_data(self):
        print 'requesting data....'
        self.response = self.engine.process(status_callback=self.update_progressbar, 
                                sources=self.sources,
                                targets=self.targets)
        print 'finished'

    def get_seismograms(self):
        return self.seismograms

    def set_seismograms(self, seismograms):
        self.seismograms = seismograms

    def set_misfit(self, misfits):
        self.misfits=misfits 

    def yaml_dump(self, fn='test.yaml'):
        '''
        Dump TestCase Object to yaml file.
        '''
        f = open(fn, 'w')
        f.write(self.dump())
        f.close()
    
    @staticmethod
    def yaml_2_TestCase(fn):
        '''
        Create a TestCase Object from file. 

        :param fn: (str) filename
        '''
        f = open(fn, 'r')
        tc = load_string(f.read())
        
    def update_progressbar(self, a, b):
        try:
            self.progressbar.update(a)
        except AttributeError:
            self.progressbar = progressbar.ProgressBar(maxval=b).start()
            self.progressbar.update(a)

    def dump_pile(self, fn='test_dumped_seismograms.mseed'):
        pile.make_pile(seismograms.values(), fn=fn)

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

    def get_sources_where(self, param_dict):
        '''
        param_dict is something like {latitude:10, longitude:10}

        :returns: list of sources, matching the required parameters in the param_dict.
        '''
        #assert all([float(v) in self.test_parameters[k] for k,v in param_dict.items()])
        return filter(lambda s: all(map(lambda x: abs(getattr(s, x[0])-x[1])<1e-7, \
                param_dict.items())), self.sources)

    def check_plot(self, param_dict):
        '''
        param_dict is something like {latitude:10, longitude:10}, defining the 
        area of source, that you would like to look at.
        '''
        sources = self.get_sources_where(param_dict)

        # google maps snuffling
        gs = gridspec.GridSpec(len(self.targets)/3,3)
        gs_dict= dict(zip(self.targets, gs))

        if self.misfit_setup.domain=='frequency_domain':
            gs_traces = gridspec.GridSpec(len(self.targets)/3,3)
            gs_traces_dict= dict(zip(self.targets, gs_traces))

        for source  in sources:
            for t,pr_cand in self.processed_candidates[source].items():
                ax = plt.subplot(gs_dict[t])
                pr_ref = self.processed_references[source][t]

                if self.misfit_setup.domain=='frequency_domain':
                    cand = pr_cand[0]
                    x = pr_cand[1]
                    y = num.log10(num.abs(pr_cand[2]))
                    ref = pr_ref[0]
                    x_ref = pr_ref[1]
                    y_ref = num.log10(num.abs(pr_ref[2]))

                    c_tracex = pr_cand[0].get_xdata()
                    c_tracey = pr_cand[0].get_ydata()

                    r_tracex = pr_ref[0].get_xdata()
                    r_tracey = pr_ref[0].get_ydata()

                    ax_t = plt.subplot(gs_traces_dict[t])
                    ax_t.set_title(t.codes, fontsize=8)
                    ax_t.plot(c_tracex, c_tracey)
                    p = ax_t.fill_between(r_tracex,
                                        0,
                                        r_tracey,
                                        facecolor='grey',
                                        alpha=0.2)

                else:
                    x = pr_cand.get_xdata()
                    y = pr_cand.get_ydata() 
                    x_ref = pr_ref.get_xdata()
                    y_ref = pr_ref.get_ydata()
                    
                ax.set_title('.'.join(t.codes), fontsize=11)
                ax.plot(x, y, label="%sW %sN %sm"%(source.lat,
                                                   source.lon, 
                                                   source.depth))

                marker_min = self.test_markers[source][t].tmin
                marker_max = self.test_markers[source][t].tmax

                plt.annotate('p', xy=(marker_min, 0))
                plt.annotate('s', xy=(marker_max, 0))
                p = ax.fill_between(x_ref,
                                    0,
                                    y_ref,
                                    facecolor='grey',
                                    alpha=0.5)
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                plt.tick_params(axis='both', which='major', labelsize=10)
                if self.misfit_setup.domain=='frequency_domain':
                    plt.xscale('log')

        plt.subplots_adjust(left=None, 
                            bottom=None, 
                            right=None, 
                            top=None, 
                            wspace=None, 
                            hspace=0.6)

        plt.legend(loc=2, prop={'size':8})
        #plt.tight_layout()
        plt.show()

    def plot1d(self, order, fix_parameter_value):
        self.numpy_it(order=order)

        x=self.num_array[0]
        y=self.num_array[1]
        z=self.num_array[2]
        v=self.num_array[3]
        X,Y,Z = griddata_auto(x,y,v)
        index = num.where(Y==fix_parameter_value)

        plt.plot(X, Z[index[0][0]])
        plt.show()

    def contourf(self, fix_parameter):
        '''
        :param fix_parameter: dict like {'lat':10}

        parameters are sorted beforehand. This also defines the x and y axis.
        (By alphabetical order)
        '''
        #sort parameters with fix_parameter key as last item
        p = self.test_parameters.keys()
        p.sort()
        p.insert(2, p.pop(p.index(fix_parameter.keys()[0])))
        self.numpy_it(order=p)

        x=self.num_array[0]
        y=self.num_array[1]
        z=self.num_array[2]
        v=self.num_array[3]
        x=x.reshape(len(self.test_parameters['lat']),
                len(self.test_parameters['lon']))
        y=y.reshape(len(self.test_parameters['lat']),
                len(self.test_parameters['lon']))
        v=v.reshape(len(self.test_parameters['lat']),
                len(self.test_parameters['lon']))
        plt.contourf(x,y,v,20,  cmap=cm.bone_r)
        plt.xlabel(self.xkey)
        plt.ylabel(self.ykey)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('L%s misfit'%int(self.misfit_setup.norm))

        plt.show()

    def ydata_of_target(self, sources, target):
        if sources == []:
            sources = self.sources
        for source in sources:
            ssmgrm = self.seismograms[source][target]
            yield source, ssmgrm.get_xdata(), ssmgrm.get_ydata()

    @property
    def store(self):
        return self.engine.get_store(self.store_id)


if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    # load stations from file:
    #stations = model.load_stations(pjoin(selfdir, '../reference_stations_castor.txt'))

    markers = gui_util.Marker.load_markers(pjoin(selfdir,
                                                    '../reference_marker_castor.txt'))

    C = Core(markers=markers, stations=None)
