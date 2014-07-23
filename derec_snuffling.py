import numpy as num
import os
import matplotlib.pyplot as plt
import copy
from collections import defaultdict
from PyQt4 import QtGui
from scipy.signal import butter

from pyrocko.gf.meta import ConfigTypeA
from pyrocko.snuffling import Snuffling, Param, Switch, NoViewerSet, Choice
from pyrocko.gf import *
from pyrocko import *
from pyrocko.guts import load_string, String, Float, Int, Dict
from pyrocko.guts_array import Array
#from pyrocko import moment_tensor_viewer as mtv

from derec import mopad
from derec import derec_utils as du
from derec import core
from derec import optics
from derec.yaml_derec import *

pjoin = os.path.join
fail_message = 'Need to load a setup, first'
derec_home = os.environ["DEREC_HOME"]
k = 1000.
misfit = None


def set_channel_if_needed(stations, traces):
    """Set channels if a station doesn't have one set. 
    Channel names are deduced from trace nslc ids"""
    for s in stations:
        if s.get_channels()==[]:
            trs = filter(lambda x:\
                    '%s.%s.%s'%(x.nslc_id[:3])==s.nsl_string(), traces)
            c_names = [t.nslc_id[3] for t in trs]
            cs = [model.Channel(name=cn) for cn in c_names]
            s.set_channels(cs)

class Derec(Snuffling):
    '''
    *pre-filter: use high/low pass filter to filter traces before processing
    
    '''
    def setup(self):
        '''Customization of the snuffling.'''
        
        self.set_name('Depth Relocation')
        self.add_parameter(Choice('Store: ', 'store_id_choice',
                                    'doctar_mainland_20Hz',
                                    ['doctar_mainland_20Hz']))

        self.add_parameter(Param('Static marker length [s]', 'static_length', \
                3, 1, 10))
        self.add_parameter(Param('Relative Marker length [%]',\
                'marker_perc_length', 0.1, 0., 100.))
        self.add_parameter(Param('Minimal Depth [m]', 'z_min', 1*k, 0., 100*k))
        self.add_parameter(Param('Maximal Depth [m]', 'z_max', 10*k, 0., 100*k))
        self.add_parameter(Param('Number of Depths', 'num_depths', 21, 1, 50))
        self.add_parameter(Param('rel. Shift of Marker []', 'marker_shift_frac',\
                0.5, 0., 1.))
        self.add_parameter(Param('Number of Time Shifts', 'num_time_shifts',\
                80., 0., 100.))
        self.add_parameter(Param('rel. trace shift [%]', 'perc_of_shift', 20.,\
                0., 100.))
        self.add_parameter(Param('Rise Time [s]', 'rise_time', 0.1, 0.1, 2.))
        #self.add_parameter(Param('Fader [s]', 'xfade',\
        #        0.3, 0.1, 10.))
        self.add_parameter(Param('Fader [s]', 'xfrac',\
                0.333, 0.01, 1.))
        self.add_parameter(Param('Highpass [Hz]', 'highpass',\
                0.7, 0.01, 100.))
        self.add_parameter(Param('Lowpass [Hz]', 'lowpass',\
                6., 1., 100.))
        self.add_parameter(Switch('Pre-filter with Main filter', 'pre_filter', False))
        self.add_trigger('Load Default Setup', self.load_setup) 
        self.add_trigger('Generate Markers', self.generate_markers) 
        self.add_trigger('Add Store Directory', self.add_store_dir)
        self.set_live_update(False)

        self.test_case_setup = None
        self.phase_ids_start = ['p','P']
        self.engine = LocalEngine(store_superdirs=[pjoin(derec_home,'fomostos')])
        self._store_ids = []
        self.sources = []
        self.targets = []
        self._ma= []
        self.__reference_source = None
        self.reference_source = None

    def call(self):
        '''Main work routine of the snuffling.'''

        self.cleanup()

        self.generate_markers()
        viewer = self.get_viewer()

        self.mts = {}
        self.best_optics = {}
    
        if not self.reference_source:
            self.active_event, self.stations = self.get_active_event_and_stations()
            self.reference_source = DCSource.from_pyrocko_event(self.active_event)




        reference_source = du.clone(self.reference_source)

        depths = num.linspace(self.z_min, self.z_max, self.num_depths)

        sources = du.test_event_generator(self.reference_source, depths)
        traces = self.chopper_selected_traces(fallback=True)
        traces = list(traces)
        traces = du.flatten_list(traces)
        if not self.targets:
            self.targets = du.stations2targets(self.stations, \
                    self.store_id_choice)
                    #measureq='HH')
        #self.get_pile().all()
        set_channel_if_needed(self.stations, traces)

        if self.pre_filter:
            traces = map(lambda x: x.copy(), traces)
            if viewer.lowpass:
                map(lambda x: x.lowpass(4, viewer.lowpass), traces)
            if viewer.highpass:
                map(lambda x: x.highpass(4, viewer.highpass), traces)
        
        self.traces_dict = du.make_traces_dict(self.reference_source, \
                self.targets,
                traces)

        stf = [[0., self.rise_time],[0.,1.]]

        norm = 2
        #taper = trace.CosFader(xfade=self.xfade)
        taper = trace.CosFader(xfrac=self.xfrac)
    
        z, p, k = butter(2, [self.lowpass*num.pi*2, self.highpass*num.pi*2.],
                           'band',
                           analog=True,
                           output='zpk')
    
        z = map(complex, z)
        p = map(complex, p)
        k = complex(k)
    
        fresponse = trace.PoleZeroResponse(z,p,k)
        fresponse.validate()
    
        misfit_setup = trace.MisfitSetup(norm=norm, 
                                         taper=taper, 
                                         domain='time_domain',
                                         filter=fresponse)
    
        test_case_setup = TestCaseSetup(reference_source=self.reference_source,
                                        sources=sources,
                                        targets=self.targets,
                                        engine=self.engine,
                                        store_id=self.store_id_choice,
                                        misfit_setup=misfit_setup,
                                        source_time_function=stf,
                                        number_of_time_shifts=int(\
                                                self.num_time_shifts),
                                        percentage_of_shift=self.perc_of_shift,
                                        static_length=self.static_length,
                                        marker_perc_length=\
                                                self.marker_perc_length,
                                        marker_shift_frac=\
                                                self.marker_shift_frac,
                                        depths=depths)
        test_case_setup.regularize()
        test_case = core.TestCase(test_case_setup)

        test_case.set_raw_references(self.traces_dict)
        test_case.set_reference_markers(self.ref_markers_dict)
        
        test_case.process(verbose=False)
        test_case.validate()
        ob = optics.OpticBase(test_case)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ob.plot_misfits(ax=ax)
        fig =plt.figure()
        ob.stack_plot()
        plt.show()


    def add_store_dir(self):
        self.engine.store_superdirs.append( str(QtGui.QFileDialog.getExistingDirectory(None, 
                                 'Open working directory', 
                                 '~', 
                                 QtGui.QFileDialog.ShowDirsOnly)))
        self._store_ids = self.engine.get_store_ids()
        self.set_parameter_choices('store_id_choice', self._store_ids) 

    def load_setup(self, fn=None):
        if fn:
            fn = self.input_filename(caption='Select a setup')
        else: 
            fn = pjoin(derec_home, 'derec', 'test_case_setup.yaml')
        f = open(fn,'r')
        self.test_case_setup = load_string(f.read())
        f.close()
        
        self.set_parameter('static_length', \
                self.test_case_setup.static_length)
        self.set_parameter('num_time_shifts', \
                float(self.test_case_setup.number_of_time_shifts))
        self.set_parameter('num_depths', \
                len(self.test_case_setup.depths))
        self.set_parameter('z_min', min(self.test_case_setup.depths))
        self.set_parameter('z_max', max(self.test_case_setup.depths))
        self.set_parameter('rise_time', \
                self.test_case_setup.source_time_function[0][1])
        self.set_parameter('marker_perc_length', \
                self.test_case_setup.marker_perc_length)
        self.set_parameter('marker_shift_frac', \
                self.test_case_setup.marker_shift_frac)
        self.set_parameter('perc_of_shift',\
                self.test_case_setup.percentage_of_shift)
        self.phase_ids_start = self.test_case_setup.phase_ids_start
        self.engine = self.test_case_setup.engine
        self.misfit_setup = self.test_case_setup.misfit_setup
        self._store_ids.extend(self.engine.get_store_ids())
        self.set_parameter_choices('store_id_choice', self._store_ids) 

        self.__orig_setup_dict = self.test_case_setup.__dict__

    def generate_markers(self):
        self.cleanup()
        self.active_event, self.stations = self.get_active_event_and_stations()

        if not self.targets:
            self.targets = du.stations2targets(self.stations, \
                    self.store_id_choice)
                    #measureq='HH')

        traces = list(self.get_pile().iter_traces())
        set_channel_if_needed(self.stations, traces)

        if not self.reference_source:
            self.reference_source = DCSource.from_pyrocko_event(self.active_event)
        
        selected_markers = self.get_viewer().selected_markers()
        selected_markers = du.make_markers_dict(self.reference_source, 
                                                self.targets,
                                                selected_markers)

        self.ref_markers_dict = du.chop_ranges(self.reference_source,
                       self.targets,
                       self.engine.get_store(self.store_id_choice),
                       self.phase_ids_start,
                       return_cache=False,
                       picked_phases=selected_markers,
                       perc=self.marker_perc_length,
                       static_length=self.static_length,
                       t_shift_frac=self.marker_shift_frac,
                       use_cake=True)
                       #channel_prefix='*')

        self._ma = self.ref_markers_dict.values()[0].values()
        
        self.add_markers(self._ma)
    
def __snufflings__():
    '''Returns a list of snufflings to be exported by this module.'''
    
    return [ Derec() ]

