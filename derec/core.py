from yaml_derec import *
from pyrocko.gf import *
from pyrocko import model, gui_util, trace, moment_tensor, io
from collections import defaultdict
from scipy.signal import butter
from pyrocko.guts import Object, load_string 
from pyrocko.guts_array import Array

import copy
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

class InvalidArguments(Exception):
    '''Is raised if an argument of the signature does not fulfull requirements. '''
    pass

def equal_attributes(o1, o2):
    '''
    Return true if two objects are equal as for their attributes. 
    '''
    return o1.__dict__ == o2.__dict__

def trace_energy(t):
    return num.sum(t.get_ydata()**2)*t.deltat

def get_gaussian_noise(nsamples, noise_scale=1.):
    t = num.linspace(0,2*num.pi,nsamples)
    n = num.cos(2*num.pi*num.random.uniform(-1,1, len(t)))*num.random.normal(
            -noise_scale,noise_scale,len(t))

    n-=n.mean()
    return n

def snr_processed(ts, tn, window, setup, pre_highpass):
    ts = ts.copy()
    tn = tn.copy()
    if window:
        tn.chop(tmin=window.tmin, tmax=window.tmax, inplace=True)
        ts.chop(tmin=window.tmin, tmax=window.tmax, inplace=True)
    ad, bd, sig_t, noi_t = ts.misfit(tn,setup, debug=True, nocache=True)
    intgr_r = trace_energy(sig_t)
    intgr_n = trace_energy(noi_t)

    SNR = todb(intgr_r/intgr_n) 
    return SNR

def snr(ts, tn, window=None, taper=None):
    ts = ts.copy()
    tn = tn.copy()
    if window:
        tn.chop(tmin=window.tmin, tmax=window.tmax, inplace=True)
        ts.chop(tmin=window.tmin, tmax=window.tmax, inplace=True)
    if taper:
        tn.taper(taper, inplace=True)
        ts.taper(taper, inplace=True)
    intgr_r = trace_energy(ts)
    intgr_n = trace_energy(tn)

    SNR = todb(intgr_r/intgr_n) 
    return SNR

def todb(val):
    return 10*num.log10(val)

def make_reference_trace(source, targets, engine, source_time_function=None,
        return_snr=False, noise_type=None, **kwargs):
    if not isinstance(source, list):
        source = [source]

    SNR = None
    sigma = None
    snr_processed = None
    sigma_pr = None
    response = engine.process(
            sources=source,
            targets=targets)
    ref_seismos = du.response_to_dict(response)
    
    if source_time_function:
        ref_seismos = du.apply_stf(ref_seismos, source_time_function)
    if noise_type is not None:
        if kwargs.get('noise', False) and noise_type=='natural':
            SNR, sigma, snr_processed, sigma_pr= natural_noise_adder(traces=ref_seismos, 
                    return_snr=return_snr, **kwargs)

        elif noise_type=='gaussian':
            SNR, sigma, snr_processed, sigma_pr= gaussian_noise_adder(traces=ref_seismos, 
                    return_snr=return_snr, **kwargs)

    if return_snr:
        return ref_seismos, (SNR, sigma), (snr_processed, sigma_pr)

    else:
        return ref_seismos

    
def gaussian_noise_adder(traces=None, noise_scale=1., setup=None, chop_ranges=None, 
        taper=None, return_snr=False, pre_highpass=None, **kwargs):
    snrs = []
    snrs_prcessed = []
    for source, targets_tr in traces.items():
        for target, tr in targets_tr.items():
            nsamples = len(tr.get_ydata())
            noise = get_gaussian_noise(nsamples, noise_scale)

            if return_snr:
                window = chop_ranges[source][target]
                rtrace = tr.copy(data=True)
                ntrace = tr.copy(data=False)
                ntrace.set_ydata(noise)
                db_ratio = snr(rtrace, ntrace, window, taper)
                if setup and pre_highpass:
                    db_ratio_processed = snr_processed(rtrace, ntrace, window, setup,
                                        pre_highpass)
                else:
                    db_ratio_processed = 999 
                snrs.append(db_ratio)
                snrs_prcessed.append(db_ratio_processed)

            tr.set_ydata(tr.get_ydata()+noise)

    N_av = num.array(snrs).mean()
    sigma = num.std(snrs)
    N_av_pr = num.array(snrs_prcessed).mean()
    sigma_pr = num.std(snrs_prcessed)
    return N_av, sigma, N_av_pr, sigma_pr


def natural_noise_adder(noise, traces, noise_scale=1., t_cut=120., chop_ranges=None, taper=None,
        return_snr=False, setup=None, pre_highpass=None, **kwargs):

    noise = make_tripets(noise)
    noise_keys = noise.keys()
    noise_target_map = {}
    snrs = []
    snrs_prcessed = []
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
                deltat_max = max(tr.deltat, noise_trace.deltat)
                tr.downsample_to(deltat_max)
                noise_trace.downsample_to(deltat_max)

            tmin_index = int(t_cut/tr.deltat)
            tmax_index = len(noise_trace.ydata)-int(t_cut/noise_trace.deltat)-len(tr.get_ydata())

            random_first_index = num.random.randint(tmin_index, tmax_index)
            last_index = random_first_index+len(tr.get_ydata())
           
            noisey = noise_trace.get_ydata()[random_first_index:last_index]
            noisey -= noisey.mean()
            noisey *= noise_scale

            if return_snr:
                window = chop_ranges[source][target]
                rtrace = tr.copy(data=True)
                ntrace = tr.copy(data=False)
                ntrace.set_ydata(noisey)
                db_ratio = snr(rtrace, ntrace, window, taper)
                db_ratio_processed = snr_processed(rtrace, ntrace, window, setup,
                        pre_highpass)
                snrs.append(db_ratio)
                snrs_prcessed.append(db_ratio_processed)

            tr.set_ydata(tr.get_ydata()+noisey)

    N_av = num.array(snrs).mean()
    sigma = num.std(snrs)
    N_av_pr = num.array(snrs_prcessed).mean()
    sigma_pr = num.std(snrs_prcessed)
    return N_av, sigma, N_av_pr, sigma_pr


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

        self.phase_cache = None
        self.blacklist = ()
        self.outliers = defaultdict()

        self.scale_minmax = False

        self.picked = None
        self.pre_highpass = None
        self.reduce_half_rise = False
        self.individual_scaling = False

        self.scaling_factors = num.linspace(0.1, 2.1, 21)
        
    def request_data(self, verbose=False):
        if verbose:
            print 'requesting data....'
            pb = self.update_progressbar
        else:
            pb = None
        try:
            self.response = self.engine.process(status_callback=pb, 
                                sources=self.sources,
                                targets=self.targets)
        except meta.OutOfBounds:
            print 'index out of bounds '
            raise

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

    def make_t_shifts(self, trac, num_samples, perc=0., seconds=0.):
        """
        :param trac: pyrocko.trace.Trace
        :param num_samples: number of time shifts
        :param perc: percentage of trace length to be shifted 
        :return: numpy array. 
        """
        if not 0. in [perc, seconds]:
            raise InvalidArguments('generating time shifted candidates requires'+\
                    'that only one of perc and seconds is given.')

        if perc:
            t_shift_max = (trac.tmax - trac.tmin) / 100. * perc
        elif seconds:
            t_shift_max = seconds

        return num.linspace(-t_shift_max/2., t_shift_max/2, num_samples)

    def make_shifted_candidates(self, source, target):
        """
        returns shifted candidates.
        """
        shifted_candidates = []
        cand = self.candidates[source][target]

        if self.test_case_setup.number_of_time_shifts==0:
            return [cand.copy()]

        t_shifts = self.make_t_shifts(cand,
                self.test_case_setup.number_of_time_shifts, 
                self.test_case_setup.percentage_of_shift,
                self.test_case_setup.time_shift)

        shifted_candidates = [cand.copy() for i in range(len(t_shifts))]
        map(lambda t,s: t.shift(s), shifted_candidates, t_shifts)
        map(lambda t: t.snap(), shifted_candidates)
        i = 0
        while i<len(shifted_candidates)-1:
            if shifted_candidates[i].tmin==shifted_candidates[i+1].tmin and \
                    shifted_candidates[i].tmax==shifted_candidates[i+1].tmax:
                        shifted_candidates.remove(shifted_candidates[i+1])
            else:
                i+=1
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
    def lines_dict(traces_dict, reduction=None, scaling=None, force_update=False):
        """
        Create matplotlib.lines.Line2D objects from traces dicts.
        :param reduction: subtract time from each x-value:
        :return lines_dict: dict with lines
        """
        scaling = {} if not scaling else scaling

        lines_dict = defaultdict(dict)
        for source, target, tr in TestCase.iter_dict(traces_dict):
            if isinstance(tr, seismosizer.SeismosizerTrace):
                tr = tr.pyrocko_trace()
            
            if reduction or force_update:
                if isinstance(reduction, float):
                    reduction_value = reduction
                elif isinstance(reduction, dict):
                    try:
                        reduction_value = reduction[source][target].tmin
                    except:
                        try:
                            reduction_value = reduction[target].tmin
                        except:
                            raise
                else:
                    try:
                        reduction_value=reduction.get_cached_tmin(target, source)
                        reduction_value+=source.time
                    except:
                        print "no valid reduction value"
                        pass
            else:
                reduction_value = 0.

            try:
                c = scaling[source]
            except KeyError:
                c = 1
            
            lines_dict[source][target] = pltlines.Line2D(
                                            tr.get_xdata()-reduction_value,
                                            tr.get_ydata()*c)

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
    
    def process(self, verbose=False, debug=False, use_cake=False):
        self.request_data(verbose)

        if verbose: print 'individual scaling c: ', self.individual_scaling

        if self.pre_highpass:
            for s,t,tr in TestCase.iter_dict(self.raw_candidates):
                tr.highpass(*self.pre_highpass)
            for s,t,tr in TestCase.iter_dict(self.raw_references):
                tr.highpass(*self.pre_highpass)

        setup = self.test_case_setup

        if verbose: print('get chopping ranges....')
        extended_test_marker, c_pc = du.chop_ranges(self.sources,
                                              self.targets,
                                              self.store,
                                              setup.phase_ids_start,
                                              return_cache=True,
                                              perc=setup.marker_perc_length,
                                              static_length=setup.static_length,
                                              t_shift_frac=\
                                                      setup.marker_shift_frac,
                                              use_cake=True)

        if self.picked or self.phase_cache:
            if verbose: print 'align phases with picked ones'
            if self.reduce_half_rise:
                if verbose: print 'and reduce by t_rise'
                static_shift = -setup.source_time_function[0][1]
            else:
                static_shift = 0.

            if self.picked:
                alignment = du.get_phase_alignment(self.picked, c_pc.as_dict)
            elif self.phase_cache:
                alignment = du.get_phase_alignment(self.phase_cache.as_dict,\
                        c_pc.as_dict)
            du.align(alignment, extended_test_marker, static_shift=static_shift)
            du.align(alignment, self.raw_candidates, static_shift=static_shift)

        self.set_candidates_markers( extended_test_marker )

        if verbose: print('chopping ref....')
        self.references = du.chop_using_markers(
                                self.raw_references, 
                                self.reference_markers, 
                                inplace=False)

        if setup.source_time_function is not None:
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
        
        self.scaled_misfits, self.scaling = self.L2_misfit(verbose=verbose,
                                                          scaling_factors=\
                                                           self.scaling_factors)

    def L2_misfit(self, verbose=False, scaling_factors=None):
        misfits, scaling = du.L2_norm(self.processed_candidates,
                             self.processed_references,
                             scaling=scaling_factors, 
                             individual_scaling=self.individual_scaling,
                             verbose=verbose)
        return misfits, scaling 


    def best_source_misfit(self, use_scaled=True, verbose=False):

        '''
        Retrieve best misfit and associated source.

        :param use_scaled: If True use misfits scaled by amplitudes.
        '''
        if verbose:
            print 'use_scaled best source: ', use_scaled

        if use_scaled:
            misfits = self.scaled_misfits
        else:
            misfits = self.misfits

        minmf = min(misfits.values())
        for s, v in misfits.iteritems():
            if v==minmf:
                return s, minmf
            

    @staticmethod
    def iter_dict(*args, **kwargs):
        """
        Iterate over a 2D-dict, yield each value.
        """
        return du.iter_dict(*args, **kwargs)



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
        event = model.Event(load=pjoin(derec_home, 'mseeds', 'doctar','doctar_2011-11-01',
                        'doctar_2011-11-01_quakefile.dat'))
        files = pjoin(derec_home, 'mseeds', 'doctar', 'doctar_2011-11-01',
        'restituted')
        traces = []
        traces.extend(io.load(fn) for fn in glob.glob(files+'/*'))
        traces = du.flatten_list(traces)
        channels = ['HHE', 'HHN', 'HHZ']

    phase_ids_start =  ['p','P']
    
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
