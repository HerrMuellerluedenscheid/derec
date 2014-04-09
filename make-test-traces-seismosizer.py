from pyrocko import util, gui_util, model, io
from tunguska import gfdb, receiver, seismosizer, source
import sys


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
        HH1  58.500 12.5000  0
        HH2  48.500 12.5000  0
        HH3  48.500  3.5000  0
        HH4  58.500  3.5000  0
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
        db = gfdb.Gfdb('fomostos/local1/local1')

        seis = seismosizer.Seismosizer(hosts=['localhost'])
        seis.set_database(db)
        seis.set_effective_dt(db.dt)
        seis.set_local_interpolation('bilinear')
        seis.set_receivers(receivers)
        seis.set_source_location(self.olat, self.olon, self.otime)
        seis.set_source_constraints(0, 0, 0, 0, 0, -1)
        self.seis = seis

    def __del__(self):
        self.seis.close()

    def __call__(self):

        # Change strike within Snuffler with the added scroll bar.
        #strike = 0
        #dip = 90
        #rake = 0
        #moment = 7.00e20
        depth = 3000
        rise_time = 1
        scale = 1E21
        mxx = 1.*scale
        mxy = 1.*scale
        myz = 1.*scale
        mxz = 1.*scale

        #explosion source
        source_params = dict(zip(['mxx', 'myy', 'mzz', 'mxy', 'mxz',
                                  'myz', 'depth', 'rise-time'],
                                 [mxx, mxx, mxx, mxy, mxz, myz,
                                  depth, rise_time]))

        s = source.Source(sourcetype='moment_tensor', sourceparams=source_params)

        #strike dip rake
        #s = source.Source('bilateral',
        #sourceparams_str ='0 0 0 %g %g %g %g %g 0 0 0 0 1 %g' % (depth, moment, strike, dip, rake, rise_time))

        self.seis.set_source(s)
        recs = self.seis.get_receivers_snapshot(which_seismograms=('syn',),
                                                which_spectra=(),
                                                which_processing='tapered')
        
        trs = []
        for rec in recs:
            for t in rec.get_traces():
                t.shift(rise_time*0.5)
                trs.append(t)

        io.save(trs, 'mseeds/%(network)s_%(station)s_%(location)s_%(channel)s.mseed')

        # Create event:
        ref_event = model.Event(lat=self.olat,
                                lon=self.olon,
                                depth=depth,
                                time=self.otime,
                                name='Reference Event')
        synthetic_event_marker = gui_util.EventMarker(event=ref_event)
        gui_util.Marker.save_markers([synthetic_event_marker], 'reference_marker.txt')


Maker = MakeTestTraces()
try:
    Maker()
except:
    Maker.seis.close()
    raise
finally:
    Maker.seis.close()

sys.exit(0)
