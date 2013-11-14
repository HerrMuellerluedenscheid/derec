from pyrocko.gui_util import PhaseMarker
from pyrocko.snuffling import Param, Snuffling, Switch
from pyrocko import cake, util
from tunguska import gfdb, receiver, seismosizer, source
import fishsod_utils as fs


class ExtendedSnuffling(Snuffling):
    def __init__(self):
        Snuffling.__init__(self)
        self.test_traces = []

    def pre_destroy(self):
        """overloaded here"""
        self.cleanup()
        if self._tempdir is not None:
            import shutil
            shutil.rmtree(self._tempdir)  
        try:
            self.my_del()
        except AttributeError:
            pass


class FindShallowSourceDepth(ExtendedSnuffling):
    """
       Find Source Depth
    """
    def __init__(self):
        ExtendedSnuffling.__init__(self)
        self.rise_time = 1.

    def my_del(self):
        """ Terminates Seismosizer"""
        if self.seis is not None:
            self.seis.close()



    def setup(self):
        
        # Give the snuffling a name:
        self.set_name('Fishsod')
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # GIVE THE GFDB DEFAULT DIRECTORY HERE:'
        gfdb_dir = 'fomostos/local1/local1'
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
        try:
            self.db = gfdb.Gfdb(gfdb_dir)
            gfdb_max_range=self.db.nx*self.db.dx-self.db.dx
            gfdb_min_range=self.db.firstx
        except OSError:
            print 'OSError: probably kiwi-tools need to be installed'
            raise
        except SystemExit:
            self.fail('Could not find Greens Functions Database at %s'%gfdb_dir)

        self.add_parameter(Param('Mxx=Myy=Mzz [Nm]', 'mxx', 1., -1., 1.))
        self.add_parameter(Param('Mxy [Nm]', 'mxy', 1., -1., 1.))
        self.add_parameter(Param('Myz [Nm]', 'myz', 1., -1., 1.))
        self.add_parameter(Param('Mxz [Nm]', 'mxz', 1., -1., 1.))
        self.add_parameter(Param('Spreading Time', 't_spread', 3., 0., 10.))
        self.add_parameter(Param('Global Time Shift', 'global_time_shift', self.rise_time/2, -1., 1.))
        self.add_parameter(Switch('Show Test Traces', 'show_test_traces', False))
        self.set_live_update(False)

    def setup_source(self, **kwargs):
        # Composition of the source
        db = self.db
        seis = seismosizer.Seismosizer(hosts=['localhost'])
        seis.set_database(db)
        seis.set_effective_dt(db.dt)
        seis.set_local_interpolation('bilinear')
        seis.set_receivers(kwargs['receivers'])
        seis.set_source_location( kwargs['origin_lat'],
                                  kwargs['origin_lon'],
                                  kwargs['otime'])
        seis.set_source_constraints(0, 0, 0, 0, 0, -1)
        self.seis = seis
        seis = None
        scale = 1E21
        source_params = dict(zip(['mxx', 'myy', 'mzz', 'mxy', 'mxz',
                                 'myz', 'depth', 'rise-time'],
                                [self.mxx*scale, self.mxx*scale, self.mxx*scale,
                                self.mxy*scale, self.mxz*scale, self.myz*scale,
                                kwargs['source_depth'], self.rise_time]))

        s = source.Source(sourcetype='moment_tensor', sourceparams=source_params)
        return s

    def call(self):
        active_event, active_stations = self.get_active_event_and_stations()
        self.cleanup()

        self.viewer = self.get_viewer()

        probe_depths = [3000, 4000]

        receivers = []

        _model = cake.load_model()
        wanted_phases = []
        wanted_phases.extend(cake.PhaseDef.classic('p'))
        phase_marker = []
        for active_station in active_stations:
            r = receiver.Receiver(lat=active_station.lat,
                                  lon=active_station.lon,
                                  depth=active_station.depth,
                                  components='neu',
                                  name='%s.%s.%s' % (active_station.network,
                                                     active_station.station,
                                                     active_station.location))

            receivers.append(r)

            rays = _model.arrivals(distances=[active_station.dist_deg],
                                   phases=wanted_phases,
                                   zstart=active_event.depth)
            for ray in rays:
                m = PhaseMarker(nslc_ids=[(active_station.network,
                                           active_station.station,
                                           '*',
                                           '*')],
                                tmin=active_event.time+ray.t+self.global_time_shift,
                                tmax=active_event.time+ray.t+self.global_time_shift+self.t_spread,
                                kind=1,
                                event=active_event,
                                incidence_angle=ray.incidence_angle(),
                                takeoff_angle=ray.takeoff_angle(),
                                phasename=ray.given_phase().definition())
                m.set_selected(True)
                phase_marker.append(m)
                self.add_marker(m)
                self.viewer.update()
        print m
        reference_pile = self.get_pile()

        test_index = 0
        for z in probe_depths:

            test_list = []
            s = self.setup_source(receivers=receivers,
                                  origin_lat=active_event.lat,
                                  origin_lon=active_event.lon,
                                  otime=active_event.time,
                                  source_depth=z)
            self.seis.set_source(s)
            try:
                recs = self.seis.get_receivers_snapshot(which_seismograms=('syn',),
                                                        which_spectra=(),
                                                        which_processing='tapered')

            except:
                print "Could not get receivers snapshot at z=%s"%z

            #probe_event = model.Event(lat=float(active_event.lat),
            #                          lon=float(active_event.lon),
            #                          depth=z,
            #                          time=active_event.time,
            #                          name='Test Event i=%s, z=%s' % (test_index, z))

            traces_to_add = []
            for rec in recs:
                for trace in rec.get_traces():
                    trace.set_codes(station='%s-%s' % (test_index, trace.station))

                    if self.show_test_traces:
                        traces_to_add.append(trace)

                    test_list.append(trace)

            #for m in self.viewer.selected_markers():
            #    print 'selected markers: ',m
            #    print m.nslc_ids
            #    chopped_reference_pile = reference_pile.chop(load_data=True,
            #                        tmin=m.tmin,
            #                        tmax=m.tmax,
            #                        trace_selector=(lambda tr: tr.nslc_id in self.viewer.selected_markers() ))
            #    print 'done'


            chopped_reference_pile = self.chopper_selected_traces()
            print 'I chopped it ', chopped_reference_pile

            #for s in chopped_reference_pile:
            #    print s
            #TODO traces_file_objects notwendiger weise laden mit chop in 2. dim
            #chopped_test_list = [tl. for tl in test_list]
            TDMF = fs.time_domain_misfit(reference_pile=chopped_reference_pile,
                                         test_list=test_list,
                                         square=True)

            FDMF = fs.frequency_domain_misfit(reference_pile=chopped_reference_pile,
                                              test_list=test_list,
                                              square=True)
            print 'time domain misfit is %s'%TDMF
            print 'frequency domain misfit is %s'%FDMF

            test_index += 1

        if self.show_test_traces:
            self.add_traces(traces_to_add)


def __snufflings__():
    return [FindShallowSourceDepth()]
