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


logger = logging.getLogger('derec_utils')


class NoMatchingTraces(Exception):
    def __str__(self):
        return 'No matching traces found.'


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
    re = 6371000.785
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
            stations.append(model.Station(network=azi_to_location_digits(azi), 
                                          station=x_to_station_digits(radius),
                                          lat=lat, 
                                          lon=lon))
            
    return stations


def z_to_network_digits(z):
    return str(int(z)/1000)[0:2].zfill(2)


def x_to_station_digits(x):
    return str(int(x)/1000)[0:3].zfill(3)


def azi_to_location_digits(azi):
    """
    :param azi:
    :return:
    """
    return str(int(azi)).zfill(3)

def make_reference_markers_cake(source, targets, model):
    
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

def chop_ranges(sources, targets, store, phase_ids_start,  phase_ids_end, 
                    static_offset=None, t_start_shift=None, t_end_shift=None, 
                    t_shift_frac=None, **kwargs):
    '''
    Create extended phase markers as preparation for chopping.

    If no phase_end value is given, takes tmax as the time of the 
    last arriving phase of phase_start.
    :return:

    static offset soll ersetzt werden....
    '''
    if kwargs.get('fallback_phases', False):
        try:
            p_fallback = kwargs['fallback_phases']['p']
        except KeyError:
            p_fallback = None
        
        try:
            s_fallback = kwargs['fallback_phases']['s']
        except KeyError:
            s_fallback = None

    assert None in [t_shift_frac, t_start_shift]
    assert None in [t_shift_frac, t_end_shift]

    phase_marker_dict = defaultdict(dict)

    for source in sources:
        for target in targets:
            dist = orthodrome.distance_accurate50m(source, target)
            args = (source.depth, dist)

            tmin = store.t('first(%s)'%phase_ids_start, args)
            if tmin is None and p_fallback:
                print 'Using p_fallback'
                tmin = store.t(p_fallback, args)
            tmin += source.time

            if static_offset:
                tmax = tmin+static_offset
            else:
                tmax = store.t('first(%s)'%phase_ids_end, args)
                if tmax is None:
                    if s_fallback is not None:
                        tmax = store.t(s_fallback, args)
                    else:
                        raise Exception("cannot interpolate tmax. Target: \n%s."%target+\
                                            '\n Source: %s'%source) 

                tmax += source.time
            
            if t_shift_frac:
                t_shift = (tmax-tmin)*t_shift_frac

            tmin -= t_shift 
            tmax += t_shift 

            m = PhaseMarker(nslc_ids=target.codes,
                            tmin=tmin,
                            tmax=tmax,
                            kind=1,
                            event=source,
                            phasename='p-s')

            m.set_selected(True)

            phase_marker_dict[source][target] = m
    return phase_marker_dict


def chop_using_markers(traces, markers, *args, **kwargs):
    '''
    Chop a list of traces or generator of traces using a list of markers.
    :rtype : list
    '''
    chopped_test_traces = defaultdict(dict)
    
    if isinstance(traces, types.GeneratorType):
        for s, t, tr in traces:
            m = markers[s][t]
            s.regularize()
            t.regularize()
            tr.chop(tmin=m.tmin,
                     tmax=m.tmax,
                     *args,
                     **kwargs)
            chopped_test_traces[s][t] = tr

    else:
        for trs in traces:
            for source, target_list in markers.items():
                for target, marker in target_list.items(): 
                    if marker.nslc_ids==trs.nslc_id:
                        trs.chop(tmin=marker.tmin,
                                 tmax=marker.tmax,
                                 *args,
                                 **kwargs)

                        chopped_test_traces[source][target]=trs

    return chopped_test_traces


def extend_phase_markers(markers=[], scaling_factor=1, stations=None, 
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
                marker.tmax = marker.get_tmin() + (marker.get_tmin() - t_event) * 0.33 * scaling_factor
                yield marker


def sampling_rate_similar(t1, t2):
    '''
    returns True, if the difference in sampling rates is bigger than 1.0% of t1's sampling rate
    '''
    return abs(t1.deltat - t2.deltat) <= t1.deltat / 100


def plot_misfit_dict(mfdict):
    plt.figure()
    print mfdict.values()
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
    candidates = test_case.seismograms
    references = test_case.references
    assert len(references.items())==1
    total_misfit = defaultdict()

    # c stuff
    #lx_norm_c = ctypes.cdll.LoadLibrary('./liblxnorm.so')
    #lx_norm_c.lxnorm_n.restype = ctypes.c_double
    #lx_norm_c.lxnorm_m.restype = ctypes.c_double
    # eo c stuff

    mfsetup = test_case.misfit_setup
    norm = mfsetup.norm

    print('calculating misfits...')
    pbar = progressbar.ProgressBar(maxval=len(sources)).start()

    for si, source in enumerate(sources):
        pbar.update(si)
        ms = num.empty([len(targets)], dtype=float)
        ns = num.empty([len(targets)], dtype=float)
        c_data = []
        r_data = []
        
        for ti, target in enumerate(targets):
            reft = references.values()[0][target]
            # hier kann man auch candidates[source].values() benutzen. geht 
            # schneller! Dafuer muessen aber erst alle candidates umsortiert werden. 
            mf = reft.misfit(candidates=[candidates[source][target]], 
                                        setups=mfsetup)
            
            for c_d, r_d , m, n in mf:
                test_case.processed_candidates[source][target] = c_d
                test_case.processed_references[source][target] = r_d
                ms[ti] = m
                ns[ti] = n

                # wichtig!
                #if reft.ydata.shape != cand.ydata.shape:
                #    raise Exception('shapes are different: %s, %s'%\
                #            (reft.ydata.shape, cand.ydata.shape))
                #
                #uydata = cand.ydata
                #vydata = reft.ydata
                ##uydata = num.random.uniform(-1e21, 1e21, len(reft.ydata))
                ##vydata = num.random.uniform(-1e21, 1e21, len(reft.ydata))
                #
                #v_c = vydata.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                #u_c = uydata.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                #norm_c = ctypes.c_double(norm)

                #size_c = ctypes.c_int(len(reft.ydata))
                ##ms[ti] = lx_norm_c.lxnorm_m(v_c, u_c, norm_c, size_c)
                ##ns[ti] = lx_norm_c.lxnorm_n(v_c, norm_c, size_c)
                ##print 'C: ', ms[ti]/ns[ti]
                #ms[ti], ns[ti] = trace.Lx_norm(uydata, vydata, norm)
                ##ms[ti], ns[ti] = trace.Lx_norm(reft.ydata, cand.ydata, norm)
                ##print 'P: ', ms[ti]/ns[ti]


        M = num.power(num.sum(abs(num.power(ms, norm))), 1./norm)
        N = num.power(num.sum(abs(num.power(ns, norm))), 1./norm)
            
        total_misfit[source] = M/N
    pbar.update(si+1)
    pbar.finish()
    test_case.set_misfit(total_misfit)
