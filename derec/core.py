from pyrocko.gf import *
from pyrocko import model, gui_util, pile, trace, moment_tensor, io
from vtkOptics import *
from collections import defaultdict
from matplotlib import cm
from matplotlib.mlab import griddata
from gmtpy import griddata_auto
from scipy.signal import butter
from guts import *
from scipy.ndimage import zoom

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
                                              phase_ids_start, 
                                              phase_ids_end)
        
        test_case.test_markers = extended_test_marker
        test_case.ref_markers = extended_ref_marker

        print('chopping ref....')
        test_case.references = du.chop_using_markers(
                                test_case.reference_seismograms.iter_results(), 
                                extended_ref_marker, 
                                t_shift_frac=0.1)

        print('chopping cand....')
        test_case.seismograms = du.chop_using_markers(
                test_case.response.iter_results(), extended_test_marker, 
                                                        t_shift_frac=0.1)

        du.calculate_misfit(test_case)
        #test_case.yaml_dump()

        # Display results===================================================
        #test_case.plot1d(order, event.lon)
        test_case.contourf(xkey='lat', ykey='lon')

        #test_case.check_plot({'lat':ref_source.lat, 'depth':ref_source.depth})

        optics = OpticBase(test_case)
        #optics.plot_1d(fix_parameters={'lat':event.lat, 'lon':event.lon})
        optics.gmt_map(stations=True, events=True)


class TestCaseSetup(Object):
    sources = List.T(Source.T()) 
    targets = List.T(Target.T()) 
    engine = Engine.T()
    store_id = String.T()
    test_parameters = List.T(String.T())
    misfit_setup = trace.MisfitSetup.T()


class TestCase(Object):
    '''
    In one test case, up to 3 parameters can be modified
    '''

    def __init__(self, test_case_setup):
        self.targets=test_case_setup.targets
        # sollte unnoetig sein:
        self.sources =test_case_setup.sources
        self.engine = test_case_setup.engine
        self.test_parameters = test_case_setup.test_parameters 
        self.store_id = test_case_setup.store_id

        self.processed_references = defaultdict(dict)
        self.references = {}
        self.processed_candidates = defaultdict(dict)
        self.candidates= {}
        self.misfits = None
        self.misfit_setup = None
        self.t_shifts = num.arange(-1, 1, 10)
            
    def request_data(self):
        print 'requesting data....'
        self.response = self.engine.process(status_callback=self.update_progressbar, 
                                sources=self.sources,
                                targets=self.targets)
        print 'finished'

    def set_candidates(self, candidates):
        self.candidates = candidates

    def set_markers(self, markers):
        self.markers = markers

    def set_misfit(self, misfits):
        self.misfits=misfits 

    def yaml_dump(self, fn='test.yaml'):
        '''
        Dump TestCase Object to yaml file.
        '''
        f = open(fn, 'w')
        f.write(self.dump())
        f.close()

    def make_shifted_candidates(self, source, target):
        """
        Generator generating shifted candidates.
        """
        for tshift in self.t_shifts:
            # needs to be copied or not?
            shifted_candidate = self.candidates[source][target]
            shifted_candidate.tmin += thift

            yield shifted_candidate

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

    def contourf(self, fix_parameter=None, xkey=None, ykey=None):
        '''
        :param fix_parameter: dict like {'lat':10}

        parameters are sorted beforehand. This also defines the x and y axis.
        (By alphabetical order)
        '''
        p = self.test_parameters.keys()
        p.sort()

        if fix_parameter:
        #sort parameters with fix_parameter key as last item
            p.insert(2, p.pop(p.index(fix_parameter.keys()[0])))
            xkey, ykey, lastkey = p

        elif xkey and ykey:
            lastkey = str((set(p)-set([xkey, ykey])).pop())
            p = [xkey, ykey, lastkey]
        self.numpy_it(order=p)

        xraw = self.num_array[0]
        yraw = self.num_array[1]
        zraw = self.num_array[2]
        vraw = self.num_array[3]

        # TODO: Ersetzen durch test parameter keys
        x=xraw.reshape(len(self.test_parameters[xkey]),
                len(self.test_parameters[ykey]))
        y=yraw.reshape(len(self.test_parameters[xkey]),
                len(self.test_parameters[ykey]))
        v=vraw.reshape(len(self.test_parameters[xkey]),
                len(self.test_parameters[ykey]))

        x = zoom(x, 4, order=1)
        y = zoom(y, 4, order=1)
        v = zoom(v, 4, order=1)
        v = num.ma.masked_where(v>2.5, v)

        palette = cm.bone_r
        cf = plt.contourf(x,y,v, 20,  cmap=cm.bone_r)

        plt.plot(getattr(self.ref_source, xkey), getattr(self.ref_source, ykey), '*')
        plt.plot(xraw, yraw, '+', color='w', markersize=4)
        plt.xlabel(self.xkey)
        plt.ylabel(self.ykey)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('L%s misfit'%int(self.misfit_setup.norm))

        plt.show()

    def ydata_of_target(self, sources, target):
        if sources == []:
            sources = self.sources
        for source in sources:
            ssmgrm = self.candidates[source][target]
            yield source, ssmgrm.get_xdata(), ssmgrm.get_ydata()

    @property
    def store(self):
        return self.engine.get_store(self.store_id)


if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    # load stations from file:
    stations = model.load_stations(pjoin(selfdir,
                            '../reference_stations_castor_selection.txt'))

    #traces = io.load(pjoin(selfdir, '../traces/2013-10-01T03-32-45/2013-10-01*'))

    markers = gui_util.Marker.load_markers(pjoin(selfdir,
                                                '../reference_marker_castor.txt'))

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

    derec_home = os.environ["DEREC_HOME"]
    store_dirs = [derec_home + '/fomostos']

    engine = LocalEngine(store_superdirs=store_dirs)
    model = get_earthmodel_from_engine(engine, store_id) 

    #TESTSOURCES===============================================
    
    zoffset= 2000.
    ref_source = du.event2source(event, 'DC', strike=37.3, dip=30, rake=-3)

    depths=num.linspace(ref_source.depth-zoffset, ref_source.depth+zoffset, 5)
    strikes = [ref_source.strike]
    dips = [ref_source.dip]
    rakes = [ref_source.rake]
    lats = [ref_source.lat]
    lons = [ref_source.lon]

    print lats, '<- lats'
    print lons, '<- lons'
    print depths, '<- depths'
    print ref_source.lat, '<- event lat'
    print ref_source.lon, '<- event lon'
    print ref_source.depth, '<- event depth'
    print event.moment_tensor.both_strike_dip_rake()[0], '<- event S D R'
    location_test_sources = [DCSource(lat=lat,
                           lon=lon,
                           depth=depth,
                           time=event.time,
                           strike=strike,
                           dip=dip,
                           rake=rake,
                           magnitude=event.magnitude) for strike in strikes 
                                                    for dip in dips 
                                                    for rake in rakes 
                                                    for lat in lats 
                                                    for lon in lons 
                                                    for depth in depths]

    for s in location_test_sources:
        s.regularize()
    #==========================================================

    reference_seismograms = make_reference_trace(ref_source,
                                                 targets, 
                                                 engine)

    # setup the misfit setup:
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

    misfit_setup = trace.MisfitSetup(norm=norm,
                                     taper=taper,
                                     domain='frequency_domain',
                                     filter=fresponse)

    test_parameters = {'depth':depths}

    test_case_setup = TestCaseSetup(sources=location_test_sources,
                                    targets=targets,
                                    engine=engine, 
                                    store_id=store_id,
                                    test_parameters=test_parameters,
                                    misfit_setup=misfit_setup)

    test_case = TestCase(test_case_setup )

    extended_ref_marker = du.chop_ranges(test_case.ref_source, 
                                        test_case.targets, 
                                        test_case.store,
                                        phase_ids_start,
                                        phase_ids_end)

    test_case.set_marker(extended_ref_marker)

    D = Doer(test_case)
