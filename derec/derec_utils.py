from math import radians, acos, sin, cos, degrees, asin, pi
import numpy as num
import logging
import matplotlib.pyplot as plt
import urllib
import types 
import ctypes 
import progressbar
from collections import defaultdict
from pyrocko import pile, util, cake, gui_util, orthodrome, trace, model
from pyrocko.gf.seismosizer import *
from pyrocko.gui_util import PhaseMarker
from pyrocko.parimap import parimap
from scipy import interpolate
from derec import phase_cache as pc

from multiprocessing import Pool, Pipe, Process, Manager
from itertools import izip

logger = logging.getLogger('derec_utils')


def make_markers_dict(source, targets, markers):
    '''
    Create the standard 2d dict like:
        dict[source][target] = marker_tmin
    '''
    targets_markers = {}
    for m in markers:
        tar = filter(lambda x: util.match_nslcs('.'.join(x.codes[:3])+'.*',
            m.nslc_ids), targets)
        if len(tar)==0:
            continue
        tar = tar[0]
        key = (source.depth, source.distance_to(tar))
        targets_markers.update({key:m.tmin-source.time})

    return targets_markers

def make_traces_dict(source, targets, traces):
    targets_traces= {}
    for tr in traces:
        tar = filter(lambda x: x.codes==tr.nslc_id, targets)
        assert len(tar)==1
        tar = tar[0]
        targets_traces.update({tar:tr})

    return {source: targets_traces}

re = 6371000.785

class NoMatchingTraces(Exception):
    def __str__(self):
        return 'No matching traces found.'


def lat_lon_relative_shift(olat, olon, north_shift=0., east_shift=0.):
    '''
    Uses flat earth approximation. 
    :param olat, olon: origin latitude/ longitude in degrees
    :param north_shift, east_shift: horizontal shifts in meters.
    '''
    rolat = radians(olat)
    dlat = north_shift/re
    dlon = east_shift/(re*cos(pi*rolat/180.))

    if north_shift==0. and east_shift!=0.:
        return olon+dlon*180/pi

    if east_shift==0. and north_shift!=0.:
        return olat+dlat*180/pi

    else: 
        return olat+dlat*180/pi, olon+dlon*180/pi


def lat_lon_from_dist_azi(olat, olon, dist, azim):
    """
    Calculate coordinates from origin latitude and longitude, distance
    and azimuth. (limitation: flatness of spheroid not included, yet)
    :param olat: Latitude of origin, type: float
    :param olon: Longitude of origin, type: float
    :param dist: Distance from origin, type: float
    :param azim: Azimuth from origin, type: float
    :return lat, lon: new latitude and longitude 
    @rtype : float, float
    """
    olat = radians(olat)
    olon = radians(olon)
    azi = radians(azim)
    b = dist/re
    a = acos(cos(b)*cos(pi/2-olat)+sin(pi/2-olat)*sin(b)*cos(azi))
    B = asin(sin(b)*sin(azi)/sin(a))
    return 90-degrees(a), degrees(B)+degrees(olon)

def station_distribution(origin, rings, **kwargs):
    '''
    origin is a tuple with (lat, lon).
    rings is a list of lists with [radius, n stations].
    '''
    rotate = {}
    if kwargs.get('rotate', False):
        rotate = kwargs['rotate']

    stations = []
    olat, olon = origin
    for radius, N in rings:
        azis = [360.*(n+1)/N for n in range(N)]

        try:
            azis = [azi + rotate[radius] for azi in azis]
            map((lambda x: x%360), azis)
        except KeyError:
            pass

        for azi in azis:
            lat, lon = lat_lon_from_dist_azi(olat, olon, radius, azi)
            stations.append(model.Station(network=str(azi_to_location_digits(azi)), 
                                          station=str(x_to_station_digits(radius)),
                                          lat=lat, 
                                          lon=lon))
            
    return stations


def z_to_network_digits(z):
    r = str(int(z)/1000)[0:2].zfill(2)
    return r


def x_to_station_digits(x):
    r = str(int(x)/1000)[0:3].zfill(3)
    return r


def azi_to_location_digits(azi):
    """
    :param azi:
    :return:
    """
    r = str(int(azi)).zfill(3)
    return r



def chop_ranges(sources, targets, store, phase_ids_start,  phase_ids_end=None,
             picked_phases={}, t_shift_frac=0., cache=True, return_cache=False, **kwargs):
    '''
    Create extended phase markers as preparation for chopping.

    If no phase_end value is given, takes tmax as the time of the 
    last arriving phase of phase_start.
    :return:
    :param static_length: length of marker
    :param perc: Percentage of tmin-t0 to add to base.

    static offset soll ersetzt werden....
    '''
    if kwargs.get('perc', False):
        perc = kwargs['perc']
    assert not phase_ids_end==perc and None in (phase_ids_end, perc)

    phase_cache = pc.PhaseCache(tmin_phase_cache=picked_phases,
                             store=store, 
                             phase_ids_start=phase_ids_start,
                             phase_ids_end=phase_ids_end)
    print 'created new phasecache instance, ', phase_cache

    parallelize = False 
    if kwargs.get('parallelize', False):
        paralellize = True

    use_cake = False
    if kwargs.get('use_cake', False):
        use_cake = True

    model = store.config.earthmodel_1d

    if not isinstance(sources, list):
        sources = [sources]

    phase_marker_dict = defaultdict(dict)
    for source in sources:
        p_event = source.pyrocko_event()
        
        for target in targets:

            tmin, tmax = phase_cache.get_cached_arrivals(target, source, **kwargs)
            
            if t_shift_frac:

                t_start_shift = abs(tmax-tmin)*t_shift_frac
                t_end_shift = t_start_shift

                tmin -= t_start_shift
                tmax -= t_end_shift 

            m = PhaseMarker(nslc_ids=[(target.codes)],
                            tmin=tmin,
                            tmax=tmax,
                            kind=1,
                            event=p_event,
                            phasename='drc')



        phase_marker_dict[source][target] = m

    if not cache:
        del phase_cache
        print 'deleted phase_cache instance'
        #phase_cache.flush()

    if return_cache:
        return phase_marker_dict, phase_cache
    else:
        return phase_marker_dict


def chop_using_markers(traces, markers, *args, **kwargs):
    '''
    Chop a list of traces or generator of traces using a list of markers as
    stencils.
    :rtype : list
    '''
    chopped_test_traces = defaultdict(dict)
    
    for source, targets_traces in traces.iteritems():
        for target, tr in targets_traces.iteritems():
            m = markers[source][target]

            chopped_tr = tr.chop(m.tmin, m.tmax, args, kwargs)
            chopped_test_traces[source][target] = chopped_tr

    return chopped_test_traces


def extend_markers(markers=[], scaling_factor=1, stations=None, 
                        event=None, phase='', model=None, inplace=False):
    '''
    Extend phase markers to fixed length proportional to time lag between
    phase markers tmin and event tmin.

    :param phase: cake.phasedef object
    '''
    distances = map(lambda s: s.dist_deg, stations)
    if event and phase and model:
        for m in markers:
            for arrivals in model.arrivals(distances=distances,
                                           phases=phase,
                                           zstart=event.depth,
                                           refine=True):
                if not inplace:
                    m = m.copy()
                m.tmax = event.time + arrivals.t
            yield m

    else:
        for marker in markers:
            if isinstance(marker, gui_util.PhaseMarker):
                t_event = marker.get_event().time
                marker.tmax = marker.get_tmin() + (marker.get_tmin() - t_event)*\
                                                            scaling_factor
                



def sampling_rate_similar(t1, t2):
    '''
    returns True, if the difference in sampling rates is bigger than 1.0% 
    of t1's sampling rate.
    '''
    return abs(t1.deltat - t2.deltat) <= t1.deltat / 100



def filter_traces_dict(self, traces_dict, tfade, freqlimits):
    for s in traces_dict.values():
        map(lambda x: x.transfer(tfade, freqlimits), s.values())


def scale_minmax(t1, t2):
    t1max = max(abs(t1.ydata))
    t2max = max(abs(t2.ydata))
    t2.ydata *= t1max/t2max


def calculate_misfit(test_case, verbose=False):
    
    sources = test_case.sources
    targets = filter(lambda x: not util.match_nslcs('.'.join(x.codes),
                                    test_case.blacklist), test_case.targets)
    references = test_case.references
    assert len(references.items())==1
    total_misfit = defaultdict()

    mfsetup = test_case.misfit_setup
    norm = mfsetup.norm

    outlier_threshold = test_case.test_case_setup.outlier_threshold

    if verbose:
        pbar = progressbar.ProgressBar(maxval=len(sources)).start()

    for si, source in enumerate(sources):
        if verbose:
            pbar.update(si)
        ms = num.empty([len(targets)], dtype=float)
        ns = num.empty([len(targets)], dtype=float)
        c_data = []
        r_data = []
        
        for ti, target in enumerate(targets):
            reft = references.values()[0][target]
            M_tmp = 999.

            shifted_candidates = test_case.make_shifted_candidates(source,
                    target)

            for cand_i in shifted_candidates:
                if test_case.scale_minmax: 
                    scale_minmax(reft, cand_i)

                m,n,r_d,c_d = reft.misfit(candidate=cand_i, setup=mfsetup, 
                        debug=True)

                if m==None or n==None:
                    print 'm,n =None'
                    continue

                if m/n>=M_tmp:
                    continue

                elif m/n<M_tmp:
                    M_tmp = m/n
                    M = m
                    N = n
                    best_candidate = c_d
                    best_reference = r_d

            try:
                test_case.processed_candidates[source][target] = best_candidate
                test_case.processed_references[source][target] = best_reference 
            except UnboundLocalError:
                print 'M > 999 in each iteration. Check data ranges!'
                continue

            ms[ti] = M
            ns[ti] = N

        M = num.power(num.sum(num.abs(num.power(ms, norm))), 1./norm)
        N = num.power(num.sum(num.abs(num.power(ns, norm))), 1./norm)
        
        M_total = M/N

        if M_total>=outlier_threshold:
            test_case.outliers[source] = M_total
        else:
            total_misfit[source] = M_total

    if verbose:
        pbar.update(si+1)
        pbar.finish()

    test_case.set_misfit(total_misfit)

def event2source(event, source_type='MT', rel_north_shift=0., rel_east_shift=0.,
        **kwargs):
    '''
    Convert pyrockos original event into seismosizer MT source.

    MT Source magnitude not scaled?!
    returns list of sources
    '''
    rel_n_deg, rel_e_deg = lat_lon_relative_shift(event.lat, event.lon,
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

        try: 
            s,d,r = kwargs['strike'], kwargs['dip'], kwargs['rake']
        except KeyError:
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
    
    return source_event


def stations2targets(stations, store_id=None, channels=[], measureq='*'):
    '''
    Convert pyrockos original stations into seismosizer targets.

    :param measureq: measuring quantity type like HH, BH or *
    '''
    targets = []
    for s in stations:
        if not channels:
            channels = s.get_channels()
        if not channels:
            channels = 'NEZ'
            
        target = [Target(codes=(s.network,s.station,s.location, measureq+component),
                                 lat=s.lat,
                                 lon=s.lon,
                                 store_id=store_id,
                                 )for component in channels]
        targets.extend(target)
    
    map(lambda x: x.regularize(), targets)
    return targets


def response_to_dict(response_dict):
    store_dict = defaultdict(dict)
    for source, target, trac in response_dict.iter_results():
        store_dict[source][target] = trac
    return store_dict


def apply_stf(traces_dict, stf):
    """
    Apply a source time function, given in timedomain.

    Rethink the time shift.......
    """
    for s, target_traces in traces_dict.items():
        for t, trac in target_traces.items():
            dt = trac.deltat
            #assert(stf[0][1]-stf[0][0] >= 2* dt)

            x = trac.get_xdata()

            y = trac.get_ydata()

            x_stf_new = num.arange(stf[0][0], stf[0][-1], dt)
            finterp = interpolate.interp1d(stf[0], stf[1])
            y_stf_new = finterp(x_stf_new)
            new_y = num.convolve(y_stf_new, y)
            trac.set_ydata(new_y)

    return traces_dict

def add_random_noise_to_trace(trace, A):
    """
    Add noise to trace, scaled with *A*.
    """
    rand = 2*num.random.random(len(trace.ydata))-1.
    trace.set_ydata(trace.get_ydata()+rand*A)


def get_tabulated_phases(engine, store_id, ids=['p']):
    """
    Extract a list of the tabulated phases within a fomosto store. 
    :param engine: pyrocko.gf.seismsizer.Engine Object
    :param store_id: string identifying the exact store to use
    :param ids: a single or a list of phase ids (strings)
    """
    if not isinstance(ids, list):
        ids = list(ids)

    ids_store = [i.id for i in  engine.get_store(store_id).config.tabulated_phases ]
    return filter(lambda x: any([x.startswith(ph) for ph in ids]), ids_store)
    

def get_earthmodel_from_engine(engine, store_id):
    """
    :return: earthmodel instance used in store of engine.
    """
    return engine.get_store(store_id).config.earthmodel_1d

def seismosizerTrace2pyrockoTrace(seismosizer_trace):
    """
    Convert a seismosizer.SeismosizerTrace into a pyrocko.trace.Trace object.
    """
    if isinstance(seismosizer_trace, trace.Trace):
        return seismosizer_trace
    t = trace.Trace(ydata=map(num.float32,seismosizer_trace.ydata), 
                        tmin=seismosizer_trace.tmin, 
                        deltat=seismosizer_trace.deltat)
    t.nslc_id = seismosizer_trace.codes
    return t

def make_lots_of_test_events(ref_source, depths, ranges, number, **kwargs):
    _sources = test_event_generator(ref_source, depths)
    return_list = []
    for i in xrange(number):

        source_list_copy = map(clone, _sources)

        set_randomized_values(source_list_copy, ranges, **kwargs)
        return_list.append(source_list_copy)

    return return_list

        
def set_randomized_values(source_list, ranges, func='uniform'):
    """
    :param func: A randomizing function like num.random.uniform (default)
    """
    values = {}
    for k,v in ranges.items():
        sourceval = getattr(source_list[0], k)
        if func=='uniform':
            val = num.random.uniform(sourceval-v, sourceval+v)
        elif func=='normal':
            val = num.random.normal(sourceval, v)
        values.update({k: val})

    for k,v in values.items():
        for source in source_list:
            setattr(source, k, v)


def test_event_generator(ref_source, depths):
    _sources = [DCSource(lat=ref_source.lat, 
                            lon=ref_source.lon, 
                            depth=depth, 
                            time=ref_source.time, 
                            strike=ref_source.strike, 
                            dip=ref_source.dip, 
                            rake=ref_source.rake, 
                            magnitude=ref_source.magnitude) for depth in depths] 
        
    map(lambda x: x.regularize(), _sources)
    return _sources


def clone(instance, **kwargs):
    d = instance.dict()
    d.update(kwargs)
    return instance.__class__(**d)


def flatten_list(the_list):
    return [item for sublist in the_list for item in sublist]
