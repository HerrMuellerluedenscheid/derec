from pyrocko.gf.seismosizer import *
from guts import *
from guts_array import *
import urllib2
from pyrocko import cake, model, gui_util, pile
import os
import derec_utils as du

pjoin = os.path.join

class Core:
    def __init__(self, markers, stations):
        base_url = 'http://kinherd.org/gfs/seismosizer?'
        store_id = 'crust2_dd' 
        
        reference_seismograms = []
        for stat in stations:
            reference_seis_req = SeismosizerRequest(store_id=store_id,
                                source_lat=0,
                                source_lon=0,
                                source_depth=2000,
                                receiver_lat=stat.lat,
                                receiver_lon=stat.lon,
                                net_code='%s' % stat.network,
                                sta_code='%s' % stat.station,
                                mnn=1.,
                                mee=1.,
                                mdd=1.,
                                mne=0.,
                                mnd=0.,
                                med=0.)

            reference_seismograms.append(request_seismogram(reference_seis_req))

        test_depths = [1000, 2000, 3000]
        test_seismograms = []
        for i, d in enumerate(test_depths):
            for stat in stations:
                reference_seis_req = SeismosizerRequest(store_id=store_id,
                                source_lat=0,
                                source_lon=0,
                                source_depth=2000,
                                receiver_lat=stat.lat,
                                receiver_lon=stat.lon,
                                net_code='%s' % stat.network,
                                sta_code='%s' % stat.station,
                                mnn=1.,
                                mee=1.,
                                mdd=1.,
                                mne=0.,
                                mnd=0.,
                                med=0.)

                test_seismograms.append(request_seismogram(reference_seis_req))

        # Request earthmodel
        model_request = request_earthmodel(store_id)
        model = model_request.store_configs[0].earthmodel_1d

        event = [m for m in markers if isinstance(m, gui_util.EventMarker)]
        assert len(event) == 1
        self.stations = stations

        # Extend P phase markers to 210 p reflection
        latest_phase = cake.PhaseDef('pPv210p')
        
        # extend picked markers
        du.extend_phase_markers(markers, phase=latest_phase)

        # chop reference_traces:
        for s in reference_seismograms:
            for t in s:
                for m in markers:
                    if util.match_nslc(t.nslc_id, m.nslc_ids):
                        print 'nefore', t
                        t.chop(tmin=m.tmin, tmax=m.tmax)
                    
                        print t
            
if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    stations = model.load_stations(pjoin(selfdir, '../reference_stations.txt')) 
    markers = gui_util.Marker.load_markers(pjoin(selfdir, '../reference_marker.txt'))
    C = Core(markers=markers, stations=stations)
