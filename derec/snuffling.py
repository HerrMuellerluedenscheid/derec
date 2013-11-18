from pyrocko.gui_util import PhaseMarker
from pyrocko.snuffling import Param, Snuffling, Switch
from pyrocko import cake, util, model
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

        _model = cake.load_model()

        receivers = map(lambda a_s:
            receiver.Receiver(lat=a_s.lat,
                                  lon=a_s.lon,
                                  depth=a_s.depth,
                                  components='neu',
                                  name='%s.%s.%s' % (a_s.network,
                                                     a_s.station,
                                                     a_s.location))

            ,active_stations)

        phase_marker = fs.phase_ranges(_model,
                                       active_stations,
                                       active_event,
                                       self.global_time_shift,
                                       self.t_spread)

        self.add_markers(phase_marker)
        self.viewer.update()
        test_index = 0
        for z in probe_depths:

            probe_event = model.Event(lat=float(active_event.lat),
                                      lon=float(active_event.lon),
                                      depth=z,
                                      time=active_event.time,
                                      name='Test Event i=%s, z=%s' % (test_index, z))

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

            traces_to_add = []
            for rec in recs:
                for trace in rec.get_traces():
                    trace.set_codes(station='%s-%s' % (test_index, trace.station))
                    test_list.append(trace)

            if self.show_test_traces:
                self.add_traces(test_list)

            probe_phase_marker = fs.phase_ranges(_model,
                                                 active_stations,
                                                 probe_event,
                                                 self.global_time_shift,
                                                 self.t_spread,
                                                 station_pref='%s-'%test_index)


            chopped_test_list = []
            for t in fs.chop_using_markers(traces=test_list, markers=probe_phase_marker):
                chopped_test_list.append(t)

            ref_selector = lambda m: m.kind == 1
            chopped_traces_groups = self.chopper_selected_traces(marker_selector=ref_selector)
            for trs in chopped_traces_groups:
                print trs

                for tr in trs:
                    print 'asdf'
                    print tr.location
                    tr.set_network('a')
                    tr.set_station('%s-%s'%(test_index, tr.station))
                    self.add_trace(tr)

            if self.show_test_traces:
                self.add_traces(traces_to_add)
                self.add_markers(probe_phase_marker)
                self.viewer.update()

            # TODO: Evtl. unterschiedliche Samplingraten beruecksichtigen!!!
            TDMF = fs.time_domain_misfit(reference_pile=chopped_traces_groups,
                                         test_list=chopped_test_list,
                                         square=True)

            FDMF = fs.frequency_domain_misfit(reference_pile=chopped_traces_groups,
                                              test_list=chopped_test_list,
                                              square=True)
            print 'time domain misfit is %s'%TDMF
            print 'frequency domain misfit is %s'%FDMF

            test_index += 1


def __snufflings__():
    return [FindShallowSourceDepth()]
