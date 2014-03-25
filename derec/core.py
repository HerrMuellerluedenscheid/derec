from pyrocko.gf import *
from pyrocko import model, gui_util, pile, trace, moment_tensor, io
from vtkOptics import *
from collections import defaultdict
from matplotlib import cm
from matplotlib.mlab import griddata
from gmtpy import griddata_auto
from scipy.signal import butter
from scipy.ndimage import zoom
from guts import *

import matplotlib.transforms as transforms
import time
import matplotlib.mlab as mlab
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.lines as pltlines
import progressbar
import os
import derec_utils as du
import numpy as num
import copy
import pdb

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
                                              perc=0.2)
                                              #phase_ids_end)
        
        test_case.set_candidates_markers( extended_test_marker )

        print('chopping ref....')
        test_case.references = du.chop_using_markers(
                                test_case.raw_references, 
                                test_case.ref_markers, 
                                t_shift_frac=0.1)

        test_case.apply_stf(test_case.test_case_setup.source_time_function)

        print('chopping cand....')
        test_case.candidates = du.chop_using_markers(
                                test_case.raw_candidates, 
                                extended_test_marker, 
                                t_shift_frac=0.1,
                                inplace=False)

        

        du.calculate_misfit(test_case)
        #test_case.yaml_dump()

        # Display results===================================================
        #test_case.plot1d(order, event.lon)
        #test_case.contourf(xkey='lat', ykey='lon')
        
        #test_case.check_plot({'lat':ref_source.lat, 'depth':ref_source.depth})

        #optics = OpticBase(test_case)
        #optics.plot_1d(fix_parameters={'lat':event.lat, 'lon':event.lon})
        #optics.gmt_map(stations=True, events=True)


class TestCaseSetup(Object):
    reference_source = Source.T()
    sources = List.T(Source.T()) 
    targets = List.T(Target.T()) 
    engine = Engine.T()
    store_id = String.T()
    test_parameters = List.T(String.T())
    misfit_setup = trace.MisfitSetup.T()
    # would be nicer in an numpy array
    source_time_function = List.T(List.T())


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
        self.misfits = None
        self.misfit_setup = test_case_setup.misfit_setup
        self.channel_map = {'N':1, 'E':2, 'Z':3}    
            
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
        self.ref_markers = markers

    def set_candidates_markers(self, markers):
        """
        candidates markers are the ones used as stencil to chop the candidates 
        traces.
        """
        self.candidates_markers = markers

    def extend_markers(self, markers, c):
        markers = du.extend_markers(markers, scaling_factor=c)
        return markers

    def targets_nsl_of(self, targets=[]):
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
        return self.targets_nsl(self.targets)

    def set_misfit(self, misfits):
        self.misfits=misfits 

    def yaml_dump(self, fn='test.yaml'):
        '''
        Dump TestCase Object to yaml file.
        '''
        f = open(fn, 'w')
        f.write(self.dump())
        f.close()

    def make_shifted_candidates(self, source, target, t_shifts):
        """
        Generator generating shifted candidates.
        """
        shifted_candidates = []
        for tshift in t_shifts:
            # needs to be copied or not?
            shifted_candidate = self.candidates[source][target].copy()
            shifted_candidate.shift(tshift)

            shifted_candidates.append(shifted_candidate)
            
        return shifted_candidates

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
        return filter(lambda s: all(map(lambda x: abs(getattr(s, x[0])-x[1])<1e-7, \
                param_dict.items())), self.sources)

    def waveforms_plot(self):
        """
        plot waveforms.
        some ideas taken from
        http://wiki.scipy.org/Cookbook/Matplotlib/MultilinePlots
        """
        num_stations = len(self.targets)/3
        figures = [plt.figure(i, facecolor='grey') for i in range(num_stations)]
        fig_dict = dict(zip(self.targets_nsl, figures))

        axes_dict = {}
        # overlapping axes objects <- google
        for target in self.targets:
            fig = fig_dict[target.codes[:3]]
            axes_dict[target] = fig.add_subplot(1,3,self.channel_map[target.codes[3]])

        processed_lines = self.lines_dict(self.processed_candidates)

        for source in self.sources:
            for target in self.targets:
                axes_dict[target].add_line(processed_lines[source][target])

        px=340.
        y_shifts = dict(zip(self.sources, num.linspace(-px/2./72., px/2./72., 
            len(self.sources))))
        
        # vertically distribute graphs
        for s in self.sources:
            for t in self.targets:
                line = processed_lines[s][t]
                fig = fig_dict[t.codes[:3]]
                trans = transforms.ScaledTranslation(0, y_shifts[s], 
                                                     fig.dpi_scale_trans)

                line.set_transform(line.get_transform()+trans)

        for ax in axes_dict.values():
            ax.axes.get_yaxis().set_visible(False)
            ax.autoscale()

    def compare_plot(self, sources=[], targets=[], traces_dicts=[], focus_first=False):
        """

        """
        def focus_data(ref_x, foc_x, foc_y):
            """
            return two lists with x- and y-data of x-range determined by *ref_x*
            """
            _foc_x = []
            _foc_y = []
            for ix,x in enumerate(ref_x):
                if foc_x[0]<=x<=foc_x[-1]:
                    _foc_y.append(foc_y[ix])
                    _foc_x.append(x)
            return _foc_x, _foc_y

        sources = self.sources if not sources else sources
        targets = self.targets if not targets else targets

        targets_nsl = self.targets_nsl_of(targets) 

        axes_dict = defaultdict(dict) 
        figs = {}

        for tnsl in targets_nsl:
            fig, axs = plt.subplots(3,1, sharex=True)
            fig.suptitle('station=%s'%'.'.join(tnsl))
            figs[tnsl] = fig

        for s in sources:
            for t in targets:
                fig = figs[t.codes[:3]]
                ax = fig.get_axes()[self.channel_map[t.codes[3]]-1]

                X = [] 
                Y = []
                for trs in traces_dicts:
                    tr = trs[s][t]
                    X.append(tr.get_xdata())
                    Y.append(tr.get_ydata())

                for i in range(len(X)-1):
                    X[i+1], Y[i+1] = focus_data(X[0], X[i+1], Y[i+1])

                [ax.plot(X[i], Y[i]) for i in range(len(X))]
                ax.set_title('%s'%t.codes[3])
                ax.autoscale()


    def lines_dict(self, traces_dict):
        """
        Create matplotlib.lines.Line2D objects from traces dicts.
        :return lines_dict: dict with lines
        """
        lines_dict = defaultdict(dict)

        for source, target_traces in traces_dict.items():
            for target, trac in target_traces.items():
                x=trac.get_xdata()
                y=trac.get_ydata()
                lines_dict[source][target] = pltlines.Line2D(x,y)

        return lines_dict 

    def stack_plot(self, param_dict=None):
        '''
        param_dict is something like {latitude:10, longitude:10}, defining the 
        area of source, that you would like to look at.
        '''
        if param_dict:
            sources = self.get_sources_where(param_dict)

        else:
            sources = self.sources

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

                marker_min = self.candidates_markers[source][t].tmin
                marker_max = self.candidates_markers[source][t].tmax

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

        plt.plot(getattr(self.reference_source, xkey), getattr(self.reference_source, ykey), '*')
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

    phase_ids_start = 'p|P'
    phase_ids_end = 's|S'

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
    
    zoffset= 0.
    ref_source = du.event2source(event, 'DC', strike=37.3, dip=30, rake=-3)

    depths=[1800, 2000, 2200]
    #depths=num.linspace(ref_source.depth-zoffset, ref_source.depth+zoffset, 1)

    print depths, '<- depths'

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
                                    source_time_function=stf)

    test_case = TestCase( test_case_setup )

    test_case.set_raw_references(reference_seismograms)

    # considering that these traces are 'real' traces. Thus, 
    # stf needs to be applied to raw traces.
    test_case.raw_references = du.apply_stf(test_case.raw_references, 
                            test_case_setup.source_time_function)

    extended_ref_marker = du.chop_ranges(ref_source, 
                                        targets, 
                                        test_case.store,
                                        phase_ids_start,
                                        perc=0.2)
                                        #phase_ids_end)

            

    test_case.set_reference_markers(extended_ref_marker)

    D = Doer(test_case)

    # Plot Distance vs. start/end of chopping range.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for source, target_marker in extended_ref_marker.items():
        for target, mark in target_marker.items():
            trac = test_case.raw_references[source][target]
            dist = source.distance_to(target)
            ax.plot(dist, mark.tmin,'+')
            ax.plot(dist, mark.tmax,'o')

    # Plot z-component for each station for each source in a different figure
    for source in test_case.sources:
        i = 0
        fig, axs = plt.subplots(len(test_case.targets)/3, sharex=True)
        sorted_targets = sorted(test_case.targets, key=lambda tr: tr.distance_to(source))
        # only use vertical components. But sort by distance, first :
        for target in [t for t in sorted_targets if t.codes[3]=='Z']:
                m = test_case.candidates_markers[source][target]
                c = test_case.raw_candidates[source][target]
                axs[i].plot(c.get_xdata(), c.get_ydata())
                axs[i].axvline(m.tmin, label='P')
                axs[i].axvline(m.tmax, label='P')

                plt.text(1, 1, '.'.join(target.codes[:3]),
                                horizontalalignment='right',
                                verticalalignment='top',
                                transform=axs[i].transAxes)
                i+=1

        fig.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
        fig.suptitle('Vertical (z) Components at z=%s'%source.depth)
    
    # plot jede Tiefe, jede station/ target: originale spur und bester candidate    

    #test_case.raw_candidates_lines = test_case.lines_dict(test_case.raw_candidates)
    #test_case.waveforms_plot()

    # TODO: traces_dicts liste erweitern um raw_references. Die haben andere
    # sources als keys!!!

    sources_z2200 = test_case.get_sources_where({'depth':2200})
    test_case.compare_plot( traces_dicts=[test_case.processed_references,
                        test_case.processed_candidates],
                        focus_first=False)
    plt.show()
    #test_case.stack_plot()
