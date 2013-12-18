from pyrocko.gf.seismosizer import *
from pyrocko import cake, model, gui_util, util
import os
import derec_utils as du
pjoin = os.path.join

class Core:
    def __init__(self, markers, stations):

        store_id = 'crust2_dd'

        event = [m for m in markers if isinstance(m, gui_util.EventMarker)]
        assert len(event) == 1
        event = event[0].get_event()

        reference_seismograms = {}
        for stat in stations:
            reference_seis_req = SeismosizerRequest(store_id=store_id,
                                                    source_lat=event.lat,
                                                    source_lon=event.lon,
                                                    source_depth=event.depth,
                                                    receiver_lat=stat.lat,
                                                    receiver_lon=stat.lon,
                                                    source_time=event.time,
                                                    net_code=stat.network,
                                                    sta_code=stat.station,
                                                    loc_code=stat.location,
                                                    mnn=1.,
                                                    mee=1.,
                                                    mdd=1.,
                                                    mne=0.,
                                                    mnd=0.,
                                                    med=0.)

            reference_seismograms[stat] = request_seismogram(reference_seis_req).traces

        test_depths = [2000, 3000]
        test_seismograms = {}
        for i, d in enumerate(test_depths):
            print 'depth: ', d
            for stat in stations:
                test_seis_req = SeismosizerRequest(store_id=store_id,
                                                   source_lat=event.lat,
                                                   source_lon=event.lon,
                                                   source_depth=d,
                                                   receiver_lat=stat.lat,
                                                   receiver_lon=stat.lon,
                                                   source_time=event.time,
                                                   net_code=stat.network,
                                                   sta_code=stat.station,
                                                   loc_code='%s-%s' % (i, stat.location),
                                                   mnn=1.,
                                                   mee=1.,
                                                   mdd=1.,
                                                   mne=0.,
                                                   mnd=0.,
                                                   med=0.)

                test_seismograms[i] = request_seismogram(test_seis_req).traces

        # Request earthmodel
        model_request = request_earthmodel(store_id)
        model = model_request.store_configs[0].earthmodel_1d

        self.stations = stations
        map(lambda s: s.set_event_relative_data(event), self.stations)

        # Extend P phase markers to 210 p reflection
        latest_phase = cake.PhaseDef('pPv210p')
        primary_phase = cake.PhaseDef('p')
        test_marker = du.chop_ranges(model, stations, event, primary_phase, latest_phase, test_depths)

        # extend picked markers
        extended_markers = list(du.extend_phase_markers(markers=markers, phase=latest_phase, stations=stations,
                                event=event, model=model))

        # chop
        chopped_ref_traces = {}
        for key, s in reference_seismograms.iteritems():
            chopped_ref_traces[key] = du.chop_using_markers(s, extended_markers)

        for d in depths:



if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    stations = model.load_stations(pjoin(selfdir, '../reference_stations.txt'))
    markers = gui_util.Marker.load_markers(pjoin(selfdir, '../reference_marker.txt'))
    print markers
    C = Core(markers=markers, stations=stations)
