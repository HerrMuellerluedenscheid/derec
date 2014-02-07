from math import radians, acos, sin, cos, degrees, asin, pi
import numpy as num
import logging
import matplotlib.pyplot as plt
import urllib
import types 
from collections import defaultdict

from pyrocko import pile, util, cake, gui_util, orthodrome, trace
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

def chop_ranges(test_case, phase_ids_start,  phase_ids_end, static_offset=None, t_start_shift=0, t_end_shift=0):
    '''
    Create extended phase markers as preparation for chopping.

    If no phase_end value is given, takes tmax as the time of the last arriving phase of phase_start.
    :return:

    static offset soll ersetzt werden....
    '''
    
    sources = test_case.sources
    targets = test_case.targets
    
    phase_marker_dict = defaultdict(dict)

    for source in sources:
        for target in targets:
            dist = orthodrome.distance_accurate50m(source, target)
            args = (source.depth, dist)
            tmin = test_case.store.t('first(%s)'%phase_ids_start, args)+source.time+t_start_shift

            if static_offset:
                tmax = tmin+static_offset
            else:
                tmax = test_case.store.t('last(%s)'%phase_ids_end, args)+source.time+t_end_shift

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


def extend_phase_markers(markers=[], scaling_factor=1, stations=None, event=None, phase='', model=None, inplace=False):
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


def calculate_misfit(test_case, mode='waveform', **kwargs):
    sources = test_case.sources
    targets = test_case.targets
    candidates = test_case.seismograms
    references = test_case.references
    assert len(references.items())==1
    total_misfit = defaultdict()

    cached_ref= {}


    mfsetup = test_case.misfit_setup

    for source in sources:
        ms = []
        ns = []
        
        for target in targets:
            reft = references.values()[0][target]
            if mode=='waveform':
                # hier kann man auch candidates[source].values() benutzen. geht 
                # schneller! Dafuer muessen aber erst alle candidates umsortiert werden. 
                mf = reft.misfit(candidates=[candidates[source][target]], 
                                            setups=mfsetup)

                for m,n in mf:
                    ms.append(m)
                    ns.append(n)

            else: 
                cand = candidates[source][target]
                cand.snap()
                reft.snap()

                try:
                    reft = cached_ref[(target, cand.tmin, cand.tmax, cand.deltat)]
                    print 'success... using a cached one '
                except KeyError:
                    max_deltat = max(cand.deltat, reft.deltat)
                    
                    if abs(reft.deltat - max_deltat) / reft.deltat > 1e-6:
                        reft.downsample_to(max_deltat, snap=True)
                    else:
                        reft.snap()

                    wanted_tmin = min(cand.tmin, reft.tmin) - max_deltat*0.5
                    wanted_tmax = max(cand.tmax, reft.tmax) + max_deltat*0.5

                    reft.extend(tmin=wanted_tmin, 
                                tmax=wanted_tmax, 
                                fillmethod='repeat')

                    
                    if kwargs.get('tfade', False):
                        tfade = kwargs[tfade]
                    else:
                        tfade = 0.0

                    if kwargs.get('freqlimits', False):
                        freqlimits = kwargs[freqlimits]
                    else:
                        freqlimits=(0.01, 0.02, 50., 100.)
                    
                    reft.transfer(tfade=tfade, 
                                  freqlimits=freqlimits, 
                                  transfer_function=mfsetup.filter)

                    if mode=='envelope':
                        reft.envelope()

                    elif mode=='positive':
                        reft.set_ydata(abs(reft.ydata))

                    cached_ref[(target, wanted_tmin, wanted_tmax, max_deltat)] =\
                                                                            reft
                if abs(reft.deltat - max_deltat) / reft.deltat > 1e-6:
                    cand.downsample_to(max_deltat, snap=True)

                cand.extend(tmin=wanted_tmin, tmax=wanted_tmax)
                cand.transfer(tfade=tfade,
                              freqlimits=freqlimits,
                              transfer_function=mfsetup.filter)
                
                if mode=='envelope':
                    cand.envelope()

                elif mode=='positive':
                    cand.set_ydata(abs(cand.ydata))

                mtmp, ntmp = trace.Lx_norm(reft.ydata, cand.ydata)
                ms.append(mtmp)
                ns.append(ntmp)
        
        ms = num.array(ms)
        ns = num.array(ns)

        norm = mfsetup.norm
        M = num.power(num.sum(num.power(ms, norm)), 1./norm)
        N = num.power(num.sum(num.power(ns, norm)), 1./norm)
            
        total_misfit[source] = M/N

    return total_misfit
