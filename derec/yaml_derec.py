import numpy as num

from guts import Object, Int, Float, Dict, List, String, Tuple, Timestamp
from guts_array import Array
from pyrocko.gf import *
from pyrocko import trace

guts_prefix = 'derec.core'

class yamlMarker(Object):
    """
    yaml version of the pyrocko.gui_util.Marker class
    """
    nslc_ids = Tuple.T(4, String.T())
    tmin = Timestamp.T()
    tmax = Timestamp.T(optional=True)
    kind = Int.T(optional=True)

class yamlTrace(Object):
    """
    To be replaced with pyrocko.gf.SeismosizerTrace
    """
    ydata = Array.T(shape=(None,), 
                    dtype=num.float32, 
                    serialize_as='base64',
                    serialize_dtype=('<f4'), 
                    optional=True)

    deltat = Float.T(optional=True)
    tmin = Float.T()
    codes = String.T()
    def __init__(self, ydata, deltat, tmin, codes):
        Object.__init__(self, ydata=ydata, deltat=deltat, tmin=tmin, codes=codes)

        def get_xdata(self):
            '''
            WARUM FUNKTIONIERT DAS NICHT???
            '''
            if self.ydata is None: raise trace.NoData()
            return self.tmin + num.arange(len(self.ydata), dtype=float64) \
                    * self.deltat

class TestCaseSetup(Object):
    tagname = String.T(default='TestCaseSetup')
    reference_source = Source.T()
    sources = List.T(Source.T()) 
    targets = List.T(Target.T()) 
    engine = Engine.T()
    store_id = String.T()
    test_parameters = List.T(String.T())
    misfit_setup = trace.MisfitSetup.T()
    # would be nicer in an numpy array
    source_time_function = List.T(List.T())
    number_of_time_shifts = Int.T()
    percentage_of_shift = Float.T()
    phase_ids_start = String.T()
    channel_map = Dict.T(String.T(), Int.T(), 
                         optional=True, 
                         default={'N':1, 'E':2, 'Z':3})


class TestCaseData(Object):
    """
    A TestCaseData object contains all waveforms, misfit results as dictionary
    and a TestCaseSetup object.
    """
    references = Dict.T(Source.T(), Dict.T(Target.T(), yamlTrace.T()), 
            optional=True)
    candidates = Dict.T(Source.T(), Dict.T(Target.T(), yamlTrace.T()), 
            optional=True)
    processed_references = Dict.T(Source.T(), Dict.T(Target.T(), yamlTrace.T()),
            optional=True)
    processed_candidates = Dict.T(Source.T(), Dict.T(Target.T(), yamlTrace.T()), 
            optional=True)

    candidates_markers = Dict.T(Source.T(), Dict.T(Target.T(), yamlMarker.T()))
    reference_markers = Dict.T(Source.T(), Dict.T(Target.T(), yamlMarker.T()))
    test_case_setup = TestCaseSetup.T(optional=True)
