from yaml_derec import *
from pyrocko.gf import *
from pyrocko import model, gui_util, trace, moment_tensor, io
from collections import defaultdict
from scipy.signal import butter
from pyrocko.guts import Object, load_string 
from pyrocko.guts_array import Array

import matplotlib.pyplot as plt
import time
import matplotlib.lines as pltlines
import progressbar
import os
import derec_utils as du
import numpy as num
import glob

pjoin = os.path.join
km = 1000.


def equal_attributes(o1, o2):
    '''
    Return true if two objects are equal as for their attributes. 
    '''
    return o1.__dict__ == o2.__dict__

def make_reference_trace(source, targets, engine, source_time_function=None, 
        noise=None):
    if not isinstance(source, list):
        source = [source]

    response = engine.process(
            sources=source,
            targets=targets)
    ref_seismos = du.response_to_dict(response)
    if source_time_function:
        ref_seismos = du.apply_stf(ref_seismos, source_time_function)
    if noise:
        noise_adder(noise, ref_seismos)

    return ref_seismos
    
def noise_adder(noise, traces, tshift=120.):
    noise = make_tripets(noise)
    noise_keys = noise.keys()
    noise_target_map = {}

    for source, targets_tr in traces.items():
        for target, tr in targets_tr.items():
            try:
                k = noise_target_map[target.codes]
            except KeyError:
                if len(noise_keys)==1:
                    k = noise_keys[0]
                    noise_keys = noise.keys()
                else:
                    k = noise_keys[num.random.randint(0, len(noise_keys))]
                try:
                    noise_keys.remove(k)
                except ValueError:
                    noise_keys = noise.keys()
                    t_shift += t_shift
                    noise_keys.remove(k)

                noise_target_map.update({target.codes: k})

            noise_trace_triplet = noise[k]
            # very unpretty:
            if 'Z' in tr.nslc_id[3]:
                noise_trace = filter(lambda x: 'Z' in x.nslc_id[3], \
                    noise_trace_triplet)
            
            elif 'E' in tr.nslc_id[3]:
                noise_trace = filter(lambda x: 'E' in x.nslc_id[3], \
                    noise_trace_triplet)
            elif 'N' in tr.nslc_id[3]:
                noise_trace = filter(lambda x: 'N' in x.nslc_id[3], \
                    noise_trace_triplet)
            assert len(noise_trace)==1
            noise_trace = noise_trace[0]
                
            if need_downsample(noise_trace, tr):
                assert noise_trace.deltat>tr.deltat
                tr.downsample_to(noise_trace.deltat)

            tmin_index = tshift*tr.deltat
            tmax_index = tmin_index+len(tr.get_ydata())

            noisey = noise_trace.get_ydata()[tmin_index:tmax_index]-\
                            noise_trace.get_ydata().mean()

            tr.set_ydata(tr.get_ydata()+noisey)


def need_downsample(t1, t2):
    return num.abs(t1.deltat-t2.deltat)>1e-4

def make_tripets(traces):
    triplets = {}
    for t in traces:
        try:
            triplets[t.nslc_id[:3]].append(t)
        except KeyError:
            triplets.update({t.nslc_id[:3]: [t]})
    return triplets


guts_prefix ='derec.yaml_derec'

class TestCase(Object):
    '''
    In one test case, up to 3 parameters can be modified
    '''
    def __init__(self, test_case_setup):
        self.set_setup(test_case_setup)

        self.raw_references = None       #(unchopped, unfiltered)
        self.processed_references = defaultdict(dict)
        self.references = {} 

        self.raw_candidates = None       #(unchopped, unfiltered)
        self.processed_candidates = defaultdict(dict)
        self.candidates= {}

    def request_data(self, verbose=False):
        if verbose:
            print 'requesting data....'
            pb = self.update_progressbar
        else:
            pb = None
        self.response = self.engine.process(status_callback=pb, 
                                sources=self.sources,
                                targets=self.targets)
        if verbose: print 'finished'
        self.set_raw_candidates(du.response_to_dict(self.response))
    
    def set_setup(self, setup=None):
        if setup:
            self.test_case_setup = setup
        self.reference_source = self.test_case_setup.reference_source
        self.targets = self.test_case_setup.targets
        self.sources = self.test_case_setup.sources
        self.engine = self.test_case_setup.engine
        self.store_id = self.test_case_setup.store_id
        self.misfit_setup = self.test_case_setup.misfit_setup
        self.channel_map = self.test_case_setup.channel_map    
        self.phase_ids_start = self.test_case_setup.phase_ids_start

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
        self.misfits = dict(misfits)

    def drop_data(self, *args):
        for k in args:
            setattr(self, k, None)

    def yaml_dump_setup(self, fn=''):
        """
        Write the setup to an individual file.
        """
        self.test_case_setup.regularize()
        f = open(fn,'w')
        dumped = self.test_case_setup.dump() 
        f.write(dumped)
        f.close()
        return dumped

    def yaml_dump(self, fn=''):
        """
        Write to yaml file.
        """

        def convert_to_yaml_dict(_dict):
            outdict = defaultdict(dict)
            for source, target_o in _dict.iteritems():
                for target, o in target_o.iteritems():
                    if isinstance(o, trace.Trace):
                        yaml_o = SeismosizerTrace.from_pyrocko_trace(o)

                    elif isinstance(o, gui_util.Marker):
                        yaml_o = yamlMarker(nslc_ids=o.nslc_ids,
                                            tmin=float(o.tmin),
                                            tmax=float(o.tmax),
                                            kind=o.kind)
                    outdict[source][target] = yaml_o

            return dict(outdict)

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

        test_case_data.misfits = dict(misfit_float_dict)

        test_case_data.reference_markers = convert_to_yaml_dict(
                self.reference_markers)

        test_case_data.candidates_markers = convert_to_yaml_dict(
                self.candidates_markers)

        test_case_data.regularize()
        test_case_data.validate()

        dumped = test_case_data.dump()
        f = open(fn, 'w')
        f.write(dumped)
        f.close()

        return dumped

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
        map(lambda t: t.snap(), shifted_candidates)
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
    def lines_dict(traces_dict, reduce=0.):
        """
        Create matplotlib.lines.Line2D objects from traces dicts.
        :param reduce: subtract time from each x-value:
        :return lines_dict: dict with lines
        """
        lines_dict = defaultdict(dict)
        reduce_value = reduce
        for source, target, tr in TestCase.iter_dict(traces_dict):
            if isinstance(tr, seismosizer.SeismosizerTrace):
                tr = tr.pyrocko_trace()

            if not isinstance(reduce, float):
                reduce_value = reduce[source][target].tmin

            lines_dict[source][target] = pltlines.Line2D(
                    tr.get_xdata()-reduce_value,
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
    
    def process(self, verbose=False, debug=False):
        self.request_data(verbose)
        setup = self.test_case_setup

        if verbose: print('chopping candidates....')
        extended_test_marker = du.chop_ranges(self.sources,
                                              self.targets,
                                              self.store,
                                              setup.phase_ids_start,
                                              perc=setup.marker_perc_length,
                                              static_length=setup.static_length,
                                              t_shift_frac=\
                                                      setup.marker_shift_frac,
                                              use_cake=True)
        
        self.set_candidates_markers( extended_test_marker )

        if verbose: print('chopping ref....')
        self.references = du.chop_using_markers(
                                self.raw_references, 
                                self.reference_markers, 
                                inplace=False)

        for s, t, tr in TestCase.iter_dict(self.references):
            tr.ydata = tr.ydata-tr.ydata.mean()

        self.apply_stf(setup.source_time_function)

        if verbose: print('chopping cand....')
        self.candidates = du.chop_using_markers(
                                self.raw_candidates, 
                                extended_test_marker, 
                                inplace=False)

        if verbose: print('calculating misfits...')

        if debug:
            trace.snuffle(self.raw_references.values()[0].values(),
                          markers=self.reference_markers.values()[0].values())

        du.calculate_misfit(self, verbose)

    def best_source_misfit(self):
        minmf = min(self.misfits.values())
        for s, v in self.misfits.iteritems():
            if v==minmf:
                return s, minmf
            

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
    store_dirs = [pjoin(derec_home, 'fomostos')]

    engine = LocalEngine(store_superdirs=store_dirs)

    test_type = 'doctar'
    add_noise = False

    if test_type == 'castor':
        store_id = 'castor'
        stations = model.load_stations(pjoin(derec_home, 'mseeds', 'castor',
                            'reference_stations_castor.txt'))
        event = model.Event(load='castor_event_2013-10-01.dat')
        channels = []
    elif test_type == 'doctar':
        #store_id = 'doctar_mainland_20Hz'
        store_id = 'castor'
        stations = model.load_stations(pjoin(derec_home, 'mseeds', 'doctar',
                        'doctar_2011-11-01', 'stations.txt'))
        event = model.Event(load=pjoin(derec_home, 'mseeds', 'doctar',
                        'doctar_2011-11-01_quakefile.dat'))
        files = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_2011-11-01',
        'restituted')
        traces = []
        traces.extend(io.load(fn) for fn in glob.glob(files+'/*'))
        traces = du.flatten_list(traces)
        channels = ['HHE', 'HHN', 'HHZ']

    phase_ids_start = '|'.join(du.get_tabulated_phases(engine,
                                                       store_id, 
                                                       ['p','P']))
    
    # load stations from file:
    # Event==================================================
    #event = filter(lambda x: isinstance(x, gui_util.EventMarker), markers)
    #assert len(event) == 1
    #event = event[0].get_event()
    #event.moment_tensor = moment_tensor.MomentTensor(strike=37.3,
    #                                                dip=30.,
    #                                                rake=-3.,
    #                                                scalar_moment=3.64e15)
    #
    #event.dump('castor_event_2013.dat')

    targets = du.stations2targets(stations, store_id, channels=channels)

    model = du.get_earthmodel_from_engine(engine, store_id) 

    #TESTSOURCES===============================================
    
    ref_source = DCSource.from_pyrocko_event(event)

    depths=[1500, 2000, 2500]

    # Das kann mit als Funktion in TestCaseSetup...
    location_test_sources = du.test_event_generator(ref_source, depths)

    map(lambda x: x.regularize(), location_test_sources)

    rise_time=1.
    
    noise = []
    if add_noise:
        noisedir = pjoin(derec_home, 'mseeds', 'iris_data', 'checked_noise')
        noise_fns = glob.glob(noisedir+'/*')
        for fn in noise_fns:
            noise.extend(io.load(fn))
        
    stf = [[0.,rise_time],[0.,1.]]
    if test_type == 'castor':
        reference_seismograms = make_reference_trace(ref_source, targets, engine,
                stf, noise=noise)

    elif test_type == 'doctar':
        reference_seismograms = du.make_traces_dict(ref_source, targets, traces)

        # not all targets have a matching trace. Thus, when iterating over targets,
        # this provokes a KeyError: 
        targets = reference_seismograms.values()[0].keys()

    # setup the misfit setup:
    norm = 2
    taper = trace.CosFader(xfrac=0.2) 
    
    z, p, k = butter(2, [0.001*num.pi*2, 2.0*num.pi*2.], 
                       'band', 
                       analog=True, 
                       output='zpk')
    
    z = map(complex, z)
    p = map(complex, p)
    k = complex(k)
    
    fresponse = trace.PoleZeroResponse(z,p,k)
    fresponse.validate()

    misfit_setup = trace.MisfitSetup(norm=norm,
                                     taper=taper,
                                     domain='time_domain',
                                     filter=fresponse)


    test_case_setup = TestCaseSetup(reference_source=ref_source,
                                    sources=location_test_sources,
                                    targets=targets,
                                    engine=engine, 
                                    store_id=store_id,
                                    misfit_setup=misfit_setup,
                                    source_time_function=stf,
                                    number_of_time_shifts=9,
                                    percentage_of_shift=10.,
                                    phase_ids_start=phase_ids_start,
                                    static_length=4.,
                                    marker_perc_length=50.,
                                    marker_shift_frac=0.2,
                                    depths=depths) 

    test_case = TestCase( test_case_setup )

    test_case.set_raw_references(reference_seismograms)

    io.save(test_case.raw_references.values()[0].values(),
            'core_traces.mseed')

    extended_ref_marker = du.chop_ranges(ref_source, 
                                        targets, 
                                        test_case.store,
                                        phase_ids_start,
                                        perc=test_case_setup.marker_perc_length,
                                        static_length=test_case_setup.static_length,
                                        t_shift_frac=test_case_setup.marker_shift_frac,
                                        use_cake=True)

    test_case.set_reference_markers(extended_ref_marker)

    test_case.process()

    plt.show()

    print 'dumping...'
    test_case.yaml_dump(fn='test_case_dump.yaml')
    print 'dumping setup'
    test_case.yaml_dump_setup(fn='test_case_setup.yaml')
