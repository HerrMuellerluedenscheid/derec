from pyrocko.gf.seismosizer import *
from pyrocko import cake, model, gui_util, util, io, pile, trace
import os
import tempfile
import derec_utils as du
import numpy as num
pjoin = os.path.join

class Core:
    def __init__(self, markers, stations):

        store_id = 'crust2_dd'

        event = [m for m in markers if isinstance(m, gui_util.EventMarker)]
        assert len(event) == 1
        event = event[0].get_event()

        reference_seismograms = []
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
                                                    mnn=0.,
                                                    mee=0.,
                                                    mdd=0.,
                                                    mne=1.,
                                                    mnd=0.5,
                                                    med=-0.5)

            reference_seismograms.extend(request_seismogram(reference_seis_req).traces)
    
        tmpdir = tempfile.mkdtemp(prefix='derec_tmp_', suffix='test')    

        io.save(reference_seismograms, filename_template=pjoin(tmpdir, 'ref.mseed'))
        test_depths = num.arange(1000,8000,500)
        test_seismograms = {}
        for i, d in enumerate(test_depths):
            seismograms = []
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
                                                   loc_code=stat.location,
                                                   mnn=0.,
                                                   mee=0.,
                                                   mdd=0.,
                                                   mne=1.,
                                                   mnd=0.,
                                                   med=0.)
                seismograms.extend(request_seismogram(test_seis_req).traces[:])
            test_seismograms[d] = seismograms

        # Check if ydata in correct depths are equal
        #for ts in test_seismograms[event.depth]:
        #    for rs in reference_seismograms:
        #        if ts.nslc_id==rs.nslc_id:
        #            assert (ts.get_ydata()==rs.get_ydata()).all()

        # Request earthmodel
        model_request = request_earthmodel(store_id)
        model = model_request.store_configs[0].earthmodel_1d

        self.stations = stations
        map(lambda s: s.set_event_relative_data(event), self.stations)

        # Extend P phase markers to 210 p reflection
        #latest_phase = cake.PhaseDef('pPv210p')
        latest_phase = cake.PhaseDef('s')
        primary_phase = cake.PhaseDef('P')
        #primary_phase = cake.PhaseDef('p')

        # markers hier ueberschreiben. Eigentlich sollen hier die gepicketn Marker verwendet werden. 

        extended_test_marker = du.chop_ranges(model, stations, event, primary_phase, latest_phase, test_depths)

        # extend picked markers

        #extended_markers = list(du.extend_phase_markers(markers=markers,
        #                                                phase=latest_phase,
        #                                                stations=stations,
        #                                                event=event, model=model))

        # chop
        chopped_ref_traces = []
        chopped_ref_traces.append(du.chop_using_markers(reference_seismograms, extended_test_marker[event.depth]))
        
        chopped_test_traces = {}
        for d in test_depths:
            chopped_test_traces[d] = du.chop_using_markers(test_seismograms[d], extended_test_marker[d]) 

        norm = 2
        
        taper = trace.CosFader(xfade=3)
        fresponse = trace.FrequencyResponse()
        setup = trace.MisfitSetup(norm=norm,
                                  taper=taper,
                                  domain='time_domain',
                                  freqlimits=(1,2,20,40),
                                  frequency_response=fresponse)
        
        total_misfit = {}
        for d,tts in test_seismograms.items():
            ms = []
            ns = []
            for rt in reference_seismograms:
                for tt in tts:
                    if rt.nslc_id==tt.nslc_id:
                        mf = rt.misfit(candidates=[tt], setups=setup)
                        for m,n in mf:
                            ms.append(m)
                            ns.append(n)
            
            ms = num.array(ms)
            ns = num.array(ns)

            M = num.sqrt(num.sum(ms**2))
            N = num.sqrt(num.sum(ns**2))
            print M, N
                
            total_misfit[d] = M/N

        print total_misfit
            
        du.plot_misfit_dict(total_misfit)

        # testweise nur element 0
        #memfile = pile.MemTracesFile(parent=None, traces=chopped_test_traces.values()[0])
        #p = pile.Pile()
        #inj = pile.Injector(p)
        #seismograms[0].snuffle()
        #inj.inject(seismograms[0])
        #from pyrocko.snuffler import snuffle
        #snuffle(memfile)

if __name__ ==  "__main__":

    selfdir = pjoin(os.getcwd(), __file__.rsplit('/', 1)[0])
    selfdir = selfdir.rsplit('/')[0]
    
    stations = model.load_stations(pjoin(selfdir, '../reference_stations.txt'))
    markers = gui_util.Marker.load_markers(pjoin(selfdir, '../reference_marker.txt'))
    C = Core(markers=markers, stations=stations)
