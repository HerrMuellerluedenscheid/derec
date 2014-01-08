from math import radians, acos, sin, cos, degrees, asin, pi
import numpy as np
import logging

from pyrocko import pile, util, cake, gui_util, orthodrome
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


def find_matching_traces(reference_pile, test_list):
    """
    Find matching pairs of traces for later Verwurstung.
    :param reference_pile:
    :param test_list:
    :return: list of sets, each containing two matching traces
    :rtype : list of sets of traces [(tr1, tr1_test) ,(tr2, tr2_test), (...),...]
    """
    trace_list = []
    num_matched = 0
    num_unmatched = 0
    for traces_group in reference_pile:
        for ref_trace in traces_group:
            for test_trace in test_list:
                if util.match_nslc('%s.%s.[0-9]-%s.%s' % (ref_trace.network,
                                                          ref_trace.station,
                                                          ref_trace.location,
                                                          ref_trace.channel),
                                   test_trace.nslc_id):
                    logger.info('Found matching traces: %s \n %s' % (ref_trace, test_trace))
                    trace_list.append([ref_trace, test_trace])
                    num_matched += 1
                    break
                else:
                    continue
            else:
                num_unmatched += 1
                logger.warning('No matching trace found for reference trace: %s' % ref_trace)
        if num_matched is 0:
            raise NoMatchingTraces()
        if num_unmatched is 0:
            logger.info('All traces found matching trace')
        else:
            logger.warning('%s of %s traces unmatched' % (num_unmatched, num_matched))

    return trace_list


def chop_longer_samples(data_set):
    """
    chop both samples to the same length.
    """
    t1 = data_set[0]
    t2 = data_set[1]
    if np.shape(t1) > np.shape(t2):
        t1 = np.resize(t1, np.shape(t2))

    elif np.shape(t1) < np.shape(t2):
        t2 = np.resize(t2, np.shape(t1))

    return [t1, t2]

def misfit_by_samples(data_set, square=False):
    """ 
    Calculate misfit from a set of two 1 dimensional arrays 

    :rtype : np.float ??? 
    :type square: bool 
    :param square: in True take square root of 
    :type data_set: set 
    :param data_set: 
    """ 
    if square: 
        exp = 2 
    else: 
        exp = 1

    data_set = chop_longer_samples(data_set)

    return np.sum(abs(pow(data_set[0], exp)-pow(data_set[1], exp)))


def frequency_domain_misfit(reference_pile, test_list, square=False, **kwargs):
    """ 

    :param reference_pile: 
    :param test_list:
    :return: 
    """ 
    traces_sets = find_matching_traces(reference_pile, test_list)
    spectra_sets = []
    for tr1, tr2 in traces_sets:
        # ignore fx-data
        spectra_sets.append(np.array((tr1.spectrum()[1],
                                      tr2.spectrum()[1])))

    return sum(map(lambda x: misfit_by_samples(x, square=square), spectra_sets))


def time_domain_misfit(reference_pile, test_list, square=False):
    """
    :param test_list:
    :param square:
    :type reference_pile: pile.Pile
    :param reference_pile:
    """
    traces_sets = find_matching_traces(reference_pile, test_list)
    downsample_if_needed(traces_sets)

    data_sets = []
    for traces_set in traces_sets:
        data_sets.append(np.array((traces_set[0].ydata, traces_set[1].ydata)))
        #map(lambda (x,y): [x.ydata, y.ydata], traces_sets)
    return sum(map(lambda x: misfit_by_samples(x, square=square), data_sets))


def chop_ranges(model, stations, event, phase_start, phase_end, depths, location_pref=''):
    '''
    Create extended phase markers as preparation for chopping.
    :param model:
    :param stations:
    :param t_spread:
    :param network_pref:
    :return:
    '''

    phase_marker_dict = {}
    phase_marker = []
    for depth in depths:
        for station in stations:
            rays = model.arrivals(distances=[station.dist_deg],
                                  phases=[phase_start, phase_end],
                                  zstart=depth,
                                  refine=True)  # how much faster is refine = false?
            rays.sort(key=lambda x: x.t)  # print dir(rays[0])
            tmin = rays[0].t
            tmax = rays[1].t
            for ray in rays:
                m = PhaseMarker(nslc_ids=[(station.network,
                                           station.station,
                                           location_pref + station.location,
                                           '*')],
                                tmin=tmin+event.time,
                                tmax=tmax+event.time,
                                kind=1,
                                event=event,
                                incidence_angle=ray.incidence_angle(),
                                takeoff_angle=ray.takeoff_angle(),
                                phasename=ray.given_phase().definition())
                m.set_selected(True)

                if location_pref:
                    m.set_kind(2)

                phase_marker.append(m)

        phase_marker_dict[depth] = phase_marker
    return phase_marker_dict


def chop_using_markers(traces, markers, *args, **kwargs):
    '''
    Chop a list of traces using a list of markers.
    :rtype : list
    '''
    chopped_test_list = []
    for marker in markers:
        for trs in traces:
            if marker.match_nslc(trs.nslc_id):
                print marker.tmin
                print marker.tmax
                print trs.tmin
                print trs.tmax
                trs.chop(tmin=marker.tmin,
                         tmax=marker.tmax,
                         *args,
                         **kwargs)

            chopped_test_list.append(trs)
    return chopped_test_list


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


def downsample_if_needed(trace_pairs):
    '''
    Downsample high frequency trace of traces pairs to the lower of both sampling rates.
    :param trace_pairs:
    :return:

    !!! resample ist problematisch, wenn die Frequenzen zu weit auseinanderliegen.
    '''

    ts = filter(lambda x: not sampling_rate_similar(x[0], x[1]), trace_pairs)
    for trace_pair in ts:
        trace_pair.sort(key=lambda x: x.deltat)
        trace_pair[0].resample(trace_pair[1].deltat)


def requests_in_gfdb_range(request, gfdb):
    '''
    Verify, that no depth in requested depths is out of range covered by the gfdb.
    '''
    if min(request) >= gfdb.firstz and max(request) <= gfdb.firstz + gfdb.nz * gfdb.dz:
        return True
