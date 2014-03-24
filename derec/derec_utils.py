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

from multiprocessing import Pool, Pipe, Process, Manager
from itertools import izip

logger = logging.getLogger('derec_utils')


re = 6371000.785

class NoMatchingTraces(Exception):
    def __str__(self):
        return 'No matching traces found.'


def lat_lon_relative_shift(olat, olon, north_shift, east_shift):
    '''
    Uses flat earth approximation. 
    :param olat, olon: origin latitude/ longitude in degrees
    :param north_shift, east_shift: horizontal shifts in meters.
    '''
    rolat = radians(olat)
    dlat = north_shift/re
    dlon = east_shift/(re*cos(pi*rolat/180.))
    
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


def cake_first_arrival(distance, depth, model, phases=None):
    """
    Get very first arrival of *phases*. 
    """
    if not phases:
        phases = ['p','P']
    
    phases = [cake.PhaseDef(ph) for ph in phases]

    tmin = min(model.arrivals([distance*cake.m2d], 
                            phases,
                            zstart=depth,
                            zstop=depth), key=lambda x: x.t).t
    return tmin


def chop_ranges(sources, targets, store, phase_ids_start,  phase_ids_end=None,
        perc=None, **kwargs):
    '''
    Create extended phase markers as preparation for chopping.

    If no phase_end value is given, takes tmax as the time of the 
    last arriving phase of phase_start.
    :return:

    static offset soll ersetzt werden....
    '''
    assert not phase_ids_end==perc and None in (phase_ids_end, perc)

    parallelize = False 
    if kwargs.get('parallelize', False):
        paralellize=True

    model = store.config.earthmodel_1d

    if not isinstance(sources, list):
        sources = [sources]

    def do_run(source, return_dict=None):
        if return_dict is None:
            return_dict = defaultdict()

        for target in targets:
            dist = source.distance_to(target)
            print dist
            args = (source.depth, dist)

            tmin = store.t('first(%s)'%phase_ids_start, args)
            if tmin==None:
                print 'tmin is None, using cake...(takes a little longer)'
                tmin = cake_first_arrival(dist, source.depth, model,
                        phases=phase_ids_start.split('|'))

            if phase_ids_end:
                tmax = store.t('first(%s)'%phase_ids_end, args)
                if tmax is None:
                    raise Exception("cannot interpolate tmax. Target: \n%s."%target+\
                                        '\n Source: %s'%source) 

            if perc:
                tmax = tmin + tmin * perc

            tmin += source.time
            tmax += source.time

            m = PhaseMarker(nslc_ids=target.codes,
                            tmin=tmin,
                            tmax=tmax,
                            kind=1,
                            event=source,
                            phasename='p-s')

            return_dict[target] = m
        return return_dict 

    phase_marker_dict = defaultdict()

    if parallelize==True:
        nworkers = 1
        #for source, tmp_dict in parimap(do_run, sources):
        manager = Manager()
        return_dict = manager.dict()
        jobs = []
        for i in range(nworkers):
            p = Process(target=do_run, args=(sources, return_dict))
            jobs.append(p)
            p.start()

        for proc in jobs:
            proc.join()
        
        #for source, tmp_dict in p.map(do_run, sources):
        #    phase_marker_dict[source] = tmp_dict
    else:
        for source in sources:
            phase_marker_dict[source] = do_run(source)

    return phase_marker_dict


def chop_using_markers(traces, markers, static_offset=None, 
                                    t_shift_frac=None, *args, **kwargs):
    '''
    Chop a list of traces or generator of traces using a list of markers.
    :rtype : list
    '''
    t_shift = 0

    chopped_test_traces = defaultdict(dict)
    
    for source, targets_traces in traces.items():
        for target, tr in targets_traces.items():
            m = markers[source][target]
            tmin = m.tmin
            tmax = m.tmax

            if static_offset:
                tmax = tmin+static_offset

            if t_shift_frac:
                t_start_shift = -(tmax-tmin)*t_shift_frac
                t_end_shift = t_start_shift

                tmin += t_start_shift
                tmax += t_end_shift 
            chopped_tr = tr.chop(tmin, tmax, args, kwargs)

            #chopped_tr.set_codes(str(source.lat), str(source.lon), str(source.depth))
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


def plot_misfit_dict(mfdict):
    plt.figure()
    plt.plot(mfdict.values(), '+')
    plt.xlabel('Depth [m]')
    plt.ylabel('Misfit []')
    plt.show()


def filter_traces_dict(self, traces_dict, tfade, freqlimits):
    for s in traces_dict.values():
        map(lambda x: x.transfer(tfade, freqlimits), s.values())


def calculate_misfit(test_case):
    
    sources = test_case.sources
    targets = test_case.targets
    references = test_case.references
    assert len(references.items())==1
    total_misfit = defaultdict()

    mfsetup = test_case.misfit_setup
    norm = mfsetup.norm

    print('calculating misfits...')
    pbar = progressbar.ProgressBar(maxval=len(sources)).start()
    import pdb
    pdb.set_trace()
    for si, source in enumerate(sources):
        pbar.update(si)
        ms = num.empty([len(targets)], dtype=float)
        ns = num.empty([len(targets)], dtype=float)
        c_data = []
        r_data = []
        
        for ti, target in enumerate(targets):
            reft = references.values()[0][target]
            #candidate = candidates[source][target]
            M_tmp = 999.

            for c_d, r_d , m, n in reft.misfit(candidates=
                        test_case.make_shifted_candidates(source, target), 
                        setups=mfsetup):
                if m==None or n==None:
                    print 'm,n =None'
                    import pdb
                    pdb.set_trace()
                    continue

                if m/n>=M_tmp:
                    continue

                elif m/n<M_tmp:
                    M_tmp = m/n
                    M = m
                    N = n
                    best_candidate = c_d
                    best_reference = r_d

            test_case.processed_candidates[source][target] = best_candidate
            test_case.processed_references[source][target] = best_reference 
            ms[ti] = M
            ns[ti] = N

        M = num.power(num.sum(num.abs(num.power(ms, norm))), 1./norm)
        N = num.power(num.sum(num.abs(num.power(ns, norm))), 1./norm)
            
        total_misfit[source] = M/N
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


def response_to_dict(response_dict):
    store_dict = defaultdict(dict)
    for source, target, trac in response_dict.iter_results():
        store_dict[source][target] = trac
    return store_dict
