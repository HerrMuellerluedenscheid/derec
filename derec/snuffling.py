from pyrocko.gui_util import PhaseMarker
from pyrocko.snuffling import Param, Snuffling, Switch
from pyrocko import io, cake, model
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
        gfdb_dir = '/data/share/u253/wegener/local2/gfdb_prem_0.05s/db'
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
        self.add_parameter(Param('STS2: Time fade [s]', 'tfade', 5., 0., 15))
        self.add_parameter(Switch('Show Test Traces', 'show_test_traces', False))
        self.add_parameter(Switch('Global Time Shift', 'add_half_rise_time', False))
        self.add_trigger('Choose GFDB', self.choose_gfdb)
        self.add_trigger('Run Taup', self.run_taup)
        self.add_trigger('Run Cake', self.run_cake)
        self.add_trigger('Save', self.savetraces)
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

        self.cleanup()

        self.viewer = self.get_viewer()
        active_event, active_stations = self.get_active_event_and_stations()
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
                                tmin=active_event.time+ray.t,
                                tmax=active_event.time+ray.t+self.t_spread,
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

                    if self.add_half_rise_time:
                        trace.tmin += self.rise_time*0.5

                    if self.show_test_traces:
                        traces_to_add.append(trace)

                    test_list.append(trace)

            chopped_reference_pile = self.chopper_selected_traces()

            #TODO traces_file_objects notwendiger weise laden mit chop in 2. dim
            chopped_test_list = [tl.chop(tmin=process_t_min, tmax=process_t_max, include_last=False) for tl in test_list]
            TDMF = fs.time_domain_misfit(reference_pile=chopped_reference_pile,
                                         test_list=chopped_test_list,
                                         square=True)
            FDMF = fs.frequency_domain_misfit(reference_pile=chopped_reference_pile,
                                              test_list=chopped_test_list,
                                              square=True)
            print 'time domain misfit is %s'%TDMF
            print 'frequency domain misfit is %s'%FDMF

            test_index += 1

        if self.show_test_traces:
            self.add_traces(traces_to_add)

        if self.auto_run_cake:
            self.run_cake()

        if self.auto_run_taup:
            self.run_taup()

    def choose_gfdb(self):
        inf = self.input_filename('Choose GFDB')
        self.db = gfdb.Gfdb(inf.rsplit('.',1)[0])
        gfdb_maxrange = self.db.nx*self.db.dx-self.db.dx
        gfdb_minrange = self.db.firstx
        self.set_parameter_range('source_depth', self.db.firstz, self.db.nz*self.db.dz-self.db.dz)
        self.set_parameter_range('dist', gfdb_minrange, gfdb_maxrange)

    def run_cake(self):
        viewer = self.get_viewer()
        for snuffling in viewer.snufflings:
            if snuffling._name is 'Cake':
                snuffling.call()
                return
        self.fail('Could not find Cake snuffling.')

    def run_taup(self):
        viewer = self.get_viewer()
        for snuffling in viewer.snufflings:
            if snuffling._name is 'TauP':
                snuffling.call()
                return
        self.fail('Could not find TauP snuffling.')
        
    def savetraces(self):
        io.save(self.test_traces, 'synthetics_depth{0}km.dist{1}km.azi{2}deg.mseed'
                .format(str(int(self.source_depth)/1000).zfill(3),
                        str(int(self.dist)/1000).zfill(3),
                        str(int(self.azi)).zfill(4)))


def __snufflings__():
    return [FindShallowSourceDepth()]
