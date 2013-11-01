from pyrocko import util, gui_util, model
from tunguska import gfdb, receiver, seismosizer, source
from numpy import array, pi, savetxt, complex, ones
import sys


class STS2:
    
    ''' Apply the STS2's transfer function which is deduced from the 
poles, zeros and gain of the transfer tunction. The Green's function database (gdfb) which is required for synthetic
seismograms and the rake of the focal mechanism can be chosen and changed within snuffler.
Two gfdbs are needed.
Three synthetic seismograms of an STS2 seismometer will be the result
'''

    def evaluate(self, freqs):
        
        # transform the frequency to angular frequency.
        w = 2j*pi*freqs

        poles = array([-3.7e-2+3.7e-2j, -3.7e-2-3.7e-2j,
                       -2.51e2, -1.31e2+4.67e2j, -1.31e2-4.67e2])
        zeros = array([0,0,0])
        k = 6.16817e7

        # Multiply factored polynomials of the transfer function's numerator
        # and denominator.
        a = ones(freqs.size,dtype=complex)*k
        for i_z in zeros:
            a *= w-i_z
        for i_p in poles:
            a /= w-i_p
        return a


def receivers_to_stations(receivers):
    stations = []
    for rec in receivers:
        stations.append(model.Station(
            network=rec.get_network(),
            station=rec.get_station(),
            location=rec.get_location(),
            lat=float(rec.lat),
            lon=float(rec.lon),
            depth=float(rec.depth)))
    return stations


class MakeTestTraces:
    def __init__(self):
        # Set up receiver configuration.
        tab = '''
        HH1  52.500  9.5000  0
        HH2  51.500  9.5000  0
        HH3  51.500  8.5000  0
        HH4  52.500  8.5000  0
        '''.strip()

        receivers = []
        for line_tab in tab.split('\n'):
            station, lat, lon, depth = line_tab.split()
            r = receiver.Receiver(lat, lon, components='neu', name='.%s.' % station)
            receivers.append(r)

        stations = receivers_to_stations(receivers)
        model.dump_stations(stations, 'reference_stations.txt')

        # Composition of the source
        self.olat, self.olon = 52.0000, 9.00000
        self.otime = util.str_to_time('1986-08-22 07:00:00')

        # The gfdb can be chosen within snuffler.
        # This refers to the 'add_parameter' method.
        db = gfdb.Gfdb('fomostos/qseis/traces')

        seis = seismosizer.Seismosizer(hosts=['localhost'])
        seis.set_database(db)
        seis.set_effective_dt(db.dt)
        seis.set_local_interpolation('bilinear')
        seis.set_receivers(receivers)
        seis.set_source_location( self.olat, self.olon, self.otime)
        seis.set_source_constraints (0, 0, 0, 0, 0, -1)
        self.seis = seis

    def __del__(self):
        self.seis.close()

    def __call__(self):

        # Change strike within Snuffler with the added scroll bar.
        strike = 0

        # Other focal mechanism parameters are constants
        dip = 90
        rake = 0
        moment = 7.00e20
        depth = 5000
        risetime = 1
        s = source.Source('bilateral',
        sourceparams_str ='0 0 0 %g %g %g %g %g 0 0 0 0 1 %g' % (depth, moment, strike, dip, rake, risetime))
        self.seis.set_source(s)
        recs = self.seis.get_receivers_snapshot( which_seismograms = ('syn',), which_spectra=(), which_processing='tapered')
        
        trs = []
        for rec in recs:
            trs.extend(rec.get_traces())
        trs.save_traces_mseed(
                filename_tmpl='mseeds/%(whichset)s_%(network)s_%(station)s_%(location)s_%(channel)s.mseed')

        # Create event:
        ref_event = model.Event(lat=self.olat,
                                lon=self.olon,
                                depth=depth,
                                time=self.otime,
                                name='Reference Event')
        synthetic_event_marker = gui_util.EventMarker(event=ref_event)
        gui_util.Marker.save_markers([synthetic_event_marker], 'reference_marker.txt')

        # Define fade in and out, band pass filter and cut off fader for the TF.
        #tfade = 8
        #freqlimit = (0.005, 0.006, 1, 1.3)
        #cut_off_fading = 5
        #ntraces = []
        #
        #for tr in trs:
        #    TF = STS2()
        #
        #    # Save synthetic trace after transfer function was applied.
        #    trace_filtered = tr.transfer(tfade, freqlimit, TF, cut_off_fading)
        #    # Set new codes to the filtered trace to make it identifiable.
        #    rename = {'e': 'BHE', 'n': 'BHN', 'u': 'BHZ'}
        #    trace_filtered.set_codes(channel=rename[trace_filtered.channel],
        #                             network='',
        #                             station='HHHA',
        #                             location='syn')
        #    ntraces.append(trace_filtered)

            ## Extract the synthetic trace's data with get_?data() and store them.
            #xval = trace_filtered.get_xdata()
            #yval = trace_filtered.get_ydata()
            #savetxt('synthetic_data_'+trace_filtered.channel,xval)

Maker = MakeTestTraces()
try:
    Maker()
except:
    Maker.seis.close()
    raise
finally:
    Maker.seis.close()

sys.exit(0)
