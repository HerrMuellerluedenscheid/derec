from numpy.core.umath import square
from pyrocko import pile, trace, util
from math import radians, acos, sin, cos, degrees, asin, pi
import numpy as np
import logging


logger = logging.getLogger('fishsod_utils')

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


def find_matching_traces(reference_pile, test_pile):
    """
    :param reference_pile:
    :param test_pile:
    :return: list of sets, each containing two matching traces
    :rtype : list
    """
    trace_list = []
    num_matched = 0
    num_unmatched = 0
    for ref_trace in reference_pile:
    #for ref_trace in reference_pile.iter_traces(load_data=True):

        # !!!MUSS!!! hier test_pile.all() sein, weil traces mit pile.add hinzugefuegt wurden!!!
        # Es sei denn, pile mit make_pile erstellt, wie im Test!!!
        # Dann muss es iter_traces(load_data=True) sein ................
        #for test_trace in test_pile.all(include_last=True):

        for test_trace in test_pile:

            #TODO: was wenn location auch mit angegeben? Wird von seismosizer auf synthetics gesetzt.
            if util.match_nslc('%s.[0-9]*%s.*.%s' % (ref_trace.network,
                                                    ref_trace.station,
                                                    ref_trace.channel),
                                                    (test_trace.nslc_id)):

                logger.info('Found matching traces: %s \n %s'%(ref_trace, test_trace))
                #TODO: evtl. copy entfernen, damit Speicher nicht volllaeuft und anders machen
                trace_list.append([ref_trace, test_trace])
                num_matched += 1
            else:                
                continue
                num_unmatched += 1
                logger.warning('No matching trace found for reference trace: %s'%ref_trace)
    if num_matched is 0:
        raise NoMatchingTraces()
    if num_unmatched is 0:
        logger.info('All traces found matching trace')
    else:
        logger.warning('%s of %s traces unmatched' % (num_unmatched, num_matched))

    return trace_list


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

    return np.sum(abs(pow(data_set[0], exp)-pow(data_set[1], exp)))


def frequency_domain_misfit(reference_pile, test_list, square=False):
    """ 

    :param reference_pile: 
    :param test_list:
    :param t_min: 
    :param t_max: 
    :return: 
    """ 
    #assert isinstance(reference_pile, pile.Pile)
    #assert isinstance(test_pile, pile.Pile)
    # convert to amp spectra:

    traces_sets = find_matching_traces(reference_pile, test_list)
    spectra_sets = []
    for tr1, tr2 in traces_sets:
        # ignore fx-data
        spectra_sets.append(np.array((tr1.spectrum()[1], tr2.spectrum()[1])))
    #map(lambda x,y:(x.spectrum(),y.spectrum()), traces_sets[0], traces_sets[1])    
    
    return sum(map(lambda x: misfit_by_samples(x, square=square), spectra_sets))


def time_domain_misfit(reference_pile, test_list, square=False):
    """

    :param test_list:
    :param square:
    :type reference_pile: pile.Pile
    :param reference_pile:
    """
    # TODO: Das ist jetzt ein tupel nach dem choppen...... Warum?
    #print type(reference_pile)
    #print type(test_pile)
    #assert isinstance(reference_pile, pile.Pile)
    #assert isinstance(test_pile, pile.Pile)

    traces_sets = find_matching_traces(reference_pile, test_list)

    #TODO: change to numpy array
    data_sets = []
    for traces_set in traces_sets:
        data_sets.append(np.array((traces_set[0].ydata, traces_set[1].ydata)))
        #map(lambda (x,y): [x.ydata, y.ydata], traces_sets)
    return sum(map(lambda x: misfit_by_samples(x, square=square), data_sets))