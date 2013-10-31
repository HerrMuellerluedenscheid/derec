from pyrocko.snuffling import Param, Snuffling, Switch
from pyrocko import util, io, model, pile
from tunguska import gfdb, receiver, seismosizer, source
from numpy import array, pi, complex, ones
import fishsod_utils as fs
import time


class STS2:
    """ Apply the STS2's transfer function which is deduced from the
        poles, zeros and gain of the transfer function. The Green's
        function database (gdfb) which is required for synthetic seismograms 
        and the rake of the focal mechanism can be chosen and changed within 
        snuffler.
        Three synthetic seismograms of an STS2 seismometer will be the result
    """
    def evaluate(self, freqs):
        
        # transform the frequency to angular frequency.
        w = 2j*pi*freqs

        poles = array([-3.7e-2+3.7e-2j, -3.7e-2-3.7e-2j,
                       -2.51e2, -1.31e2+4.67e2j, -1.31e2-4.67e2])
        #2 zeros -> displacement
        #3 zeros -> velocity
        zeros = array([0, 0, 0])
        k = 6.16817e7

        # Multiply factored polynomials of the transfer function's numerator
        # and denominator.
        a = ones(freqs.size, dtype=complex)*k
        for i_z in zeros:
            a *= w-i_z
        for i_p in poles:
            a /= w-i_p
        return a


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
        seis.set_source_constraints (0, 0, 0, 0, 0, -1)
        self.seis = seis
        seis = None

        source_params = dict(zip(['mxx', 'myy', 'mzz', 'mxy', 'mxz',
                                 'myz', 'depth', 'rise-time'],
                                [self.mxx, self.mxx, self.mxx,
                                self.mxy, self.mxz, self.myz,
                                kwargs['source_depth'], self.rise_time]))

        s = source.Source(sourcetype='moment_tensor', sourceparams=source_params)
        return s

    def setup(self):
        
        # Give the snuffling a name:
        self.set_name('Fishsod')
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # GIVE THE GFDB DEFAULT DIRECTORY HERE:'
        gfdb_dir = '/data/share/u253/wegener/local2/gfdb/gemini-iasp91-20000km/db'
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

        # Add scroll bars of the parameters that you desire to adjust.
        # 1st argument: Description that appears within the snuffling.
        # 2nd argument: Name of parameter as used in the following code.
        # 3rd-5th argument: default, start, stop.

        #self.add_parameter(Param('Source Depth [m]',
        #                         'source_depth',
        #                         self.db.firstz,
        #                         self.db.firstz,
        #                         self.db.nz*self.db.dz-self.db.dz))
        #self.add_parameter(Param('Distance [m]', 'dist', gfdb_max_range/2,
        #                         gfdb_min_range, gfdb_max_range))
        #self.add_parameter(Param('Azimuth [deg]', 'azi', 0., -180., 180.))
        self.add_parameter(Param('Mxx=Myy=Mzz [Nm]', 'mxx', 1., -1., 1.))
        self.add_parameter(Param('Mxy [Nm]', 'mxy', 0., -1., 1.))
        self.add_parameter(Param('Myz [Nm]', 'myz', 0., -1., 1.))
        self.add_parameter(Param('Mxz [Nm]', 'mxz', 0., -1., 1.))
        self.add_parameter(Param('STS2: Time fade [s]', 'tfade', 5., 0., 15))
        self.add_parameter(Switch('simulate STS2', 'simulate_STS2', False))
        self.add_parameter(Switch('Show Test Traces', 'show_test_traces', False))
        self.add_parameter(Switch('autorun Taup', 'auto_run_taup', False))
        self.add_parameter(Switch('autorun Cake', 'auto_run_cake', False))
        self.add_parameter(Switch('Add 0.5*rise time', 'add_half_rise_time', True))
        self.add_trigger('Choose GFDB', self.choose_gfdb)
        self.add_trigger('Run Taup', self.run_taup)
        self.add_trigger('Run Cake', self.run_cake)
        self.add_trigger('Save', self.savetraces)
        self.set_live_update(False)

    def call(self):

        self.cleanup()
        
        self.viewer = self.get_viewer()
        #print dir(self.viewer.get_panel_parent())

        active_event, active_stations = self.get_active_event_and_stations()

        probe_depths = [1000, 5000, 10000, 50000]

        receivers = []

        for active_station in active_stations:
            r = receiver.Receiver(lat=active_station.lat,
                                  lon=active_station.lon,
                                  depth=active_station.depth,
                                  components='neu',
                                  name='%s.%s.%s' % (active_station.network,
                                                     active_station.station,
                                                     active_station.location))

            receivers.append(r)

        reference_pile = self.get_pile()

        process_t_min = active_event.time+0.5
        process_t_max = active_event.time+20

        test_index = 0
        for z in probe_depths:

            #test_list = pile.TracesFile()
            test_list = []
            s = self.setup_source(receivers=receivers,
                                  origin_lat=active_event.lat,
                                  origin_lon=active_event.lon,
                                  otime=active_event.time,
                                  source_depth=z)
            self.seis.set_source(s)
            recs = self.seis.get_receivers_snapshot(which_seismograms=('syn',),
                                                    which_spectra=(),
                                                    which_processing='tapered')

            probe_event = model.Event(lat=float(active_event.lat),
                                      lon=float(active_event.lon),
                                      depth=z,
                                      time=active_event.time,
                                      name='Test Event i=%s, z=%s' % (test_index, z))
            #self.viewer.add_event(probe_event)

            #rename = {'e': 'BHE', 'n': 'BHN', 'u': 'BHZ'}
            for rec in recs:
                for trace in rec.get_traces():
                    trace.set_codes(station='%s-%s' % (test_index, trace.station))

                    if self.add_half_rise_time:
                        trace.tmin += self.rise_time*0.5

                    #if self.show_test_traces:
                    #    self.add_traces(trace)

                    test_list.append(trace)

            chopped_reference_pile, tracesFileObjects = reference_pile.chop(tmin=process_t_min,
                                                                            tmax=process_t_max,
                                                                            load_data=True)
            #TODO traces_file_objects notwendiger weise laden mit chop in 2. dim
            chopped_test_list = [tl.chop(tmin=process_t_min, tmax=process_t_max, include_last=True) for tl in test_list]
            TDMF = fs.time_domain_misfit(reference_pile=chopped_reference_pile,
                                         test_list=chopped_test_list,
                                         square=True)
            FDMF = fs.frequency_domain_misfit(reference_pile=chopped_reference_pile,
                                              test_list=chopped_test_list,
                                              square=True)
            print 'time domain misfit is %s'%TDMF
            print 'frequency domain misfit is %s'%FDMF

            test_index += 1

        #if self.simulate_STS2:
        #    # Define fade in and out, band pass filter and cut off fader for the TF.
        #    freq_limit = (0.005, .006, 1, 1.2)
        #    cut_off_fading = 300
        #    ntraces = []
        #
        #    for tr in self.test_traces:
        #        TF = STS2()
        #        trace_filtered = tr.transfer(self.tfade, freq_limit, TF, cut_off_fading)
        #        # Set new codes to the filtered trace to make it identifiable.
        #
        #        trace_filtered.set_codes(channel=rename[trace_filtered.channel],
        #                                 network='STS2',
        #                                 station='HH',
        #                                 location='syn')
        #    ntraces.append(trace_filtered)
        #    self.add_traces(ntraces)

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
