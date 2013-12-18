from pyrocko.snuffling import Param, Snuffling, Switch
from pyrocko import cake, model, mopad
from tunguska import gfdb, receiver, seismosizer, source
import derec_utils as du
import matplotlib.pyplot as plt


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
        self.set_name('Derec')
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # GIVE THE GFDB DEFAULT DIRECTORY HERE:'
        gfdb_dir = 'fomostos/local1/local1'
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        try:
            self.db = gfdb.Gfdb(gfdb_dir)
        except OSError:
            print 'OSError: probably kiwi-tools need to be installed'
            raise
        except SystemExit:
            self.fail('Could not find Greens Functions Database at %s' % gfdb_dir)

        self.add_parameter(Param('Mxx=Myy=Mzz [Nm]', 'mxx', 1., -1., 1.))
        self.add_parameter(Param('Mxy [Nm]', 'mxy', 1., -1., 1.))
        self.add_parameter(Param('Myz [Nm]', 'myz', 1., -1., 1.))
        self.add_parameter(Param('Mxz [Nm]', 'mxz', 1., -1., 1.))
        self.add_parameter(Param('Spreading Time', 't_spread', 3., 0., 10.))
        #self.add_parameter(Param('Global Time Shift', 'global_time_shift', self.rise_time/2, -1., 1.))
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
        seis.set_source_location(kwargs['origin_lat'],
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
        active_event, active_stations = self.get_active_event_and_stations()
        viewer = self.get_viewer()

        probe_depths = [2000,2500, 3000, 3500 ,4000,5000]
        results = {}
        TDMF = []
        FDMF = []

        if not du.requests_in_gfdb_range(probe_depths, self.db):
            self.fail('GFDB doesn\'t cover requested depth range.')

        _model = cake.load_model()

        receivers = map(lambda a_s:
                        receiver.Receiver(lat=a_s.lat,
                                          lon=a_s.lon,
                                          depth=a_s.depth,
                                          components='neu',
                                          name='%s.%s.%s' % (a_s.network,
                                                             a_s.station,
                                                             a_s.location))

            , active_stations)

        du.extend_phase_markers(viewer.markers)
        viewer.update()

        # Oder mit itertools.tee()?! mal schauen...
        chopped_reference_groups = list(self.chopper_selected_traces())

        test_index = 0
        for z in probe_depths:

            probe_event = model.Event(lat=float(active_event.lat),
                                      lon=float(active_event.lon),
                                      depth=z,
                                      time=active_event.time,
                                      name='Test Event i=%s, z=%s' % (test_index, z))

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
            #TODO: logger warning abfangen.
            except:
                print "Could not get receivers snapshot at z=%s" % z
                raise

            test_list = []
            for rec in recs:
                for trace in rec.get_traces():
                    trace.shift(self.rise_time * 0.5)
                    trace.set_codes(network='%s-%s' % (test_index, trace.network),
                                    location='sym')
                    test_list.append(trace)

            probe_phase_marker = du.chop_ranges(_model, active_stations, probe_event, self.t_spread,
                                                t_spread='%s-' % test_index)

            chopped_test_list = du.chop_using_markers(traces=test_list, markers=probe_phase_marker)

            if self.show_test_traces:

                for trs in chopped_reference_groups:
                    for tr in trs:
                        self.add_trace(tr)

                self.add_traces(chopped_test_list)
                self.add_markers(probe_phase_marker)
                viewer.update()

            # TODO: Evtl. unterschiedliche Samplingraten beruecksichtigen!!!
            TDMF.append(du.time_domain_misfit(reference_pile=chopped_reference_groups,
                                              test_list=chopped_test_list,
                                              square=True))

            FDMF.append(du.frequency_domain_misfit(reference_pile=chopped_reference_groups,
                                                   test_list=chopped_test_list,
                                                   square=True))

            test_index += 1

        results['tdmf'] = TDMF
        results['fdmf'] = FDMF

        #fig =self.figure(name='asdf')
        #fig = self.pylab(get='figure')
        plt.figure()
        plt.subplot(211)
        plt.plot(probe_depths, results['tdmf'], '-o')
        plt.title('Time Domain')
        plt.xlabel('Depth')
        plt.subplot(212)
        plt.plot(probe_depths, results['fdmf'], '-o')
        plt.title('Frequency Domain')
        plt.xlabel('Depth')

        #M = mopad.MomentTensor([1,1,1,1,1,1])
        #bb = mopad.BeachBall(M)
        #bb.ploBB({})
        plt.show()

def __snufflings__():
    return [FindShallowSourceDepth()]
