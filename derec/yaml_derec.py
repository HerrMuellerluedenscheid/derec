import numpy as num

from pyrocko.guts import Object, Int, Float, List, String, Tuple, Timestamp, Dict
from pyrocko.guts_array import Array
from pyrocko.gf import *
from pyrocko import trace

guts_prefix = 'derec.yaml_derec'

class yamlMarker(Object):
    """
    yaml version of the pyrocko.gui_util.Marker class
    """
    nslc_ids = Tuple.T(4, String.T())
    tmin = Timestamp.T()
    tmax = Timestamp.T(optional=True)
    kind = Int.T(optional=True)

    def __init__(self, nslc_ids, tmin, tmax, kind):
        Object.__init__(self, nslc_ids=nslc_ids, tmin=tmin, tmax=tmax,
                kind=kind)
        
    @staticmethod
    def from_pyrocko_marker(self, marker):
        """
        Create yamlMarker instance from pyrocko.gui_util.Marker instance.
        """
        return self.__init__(marker.nslc_ids, marker.tmin, \
                marker.tmax, marker.kind)

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
        if self.ydata is None: raise trace.NoData()
        return self.tmin + num.arange(len(self.ydata), dtype=num.float64) \
                * self.deltat

    def get_ydata(self):
        return self.ydata

class TestCaseSetup(Object):
    tagname = String.T(default='TestCaseSetup')
    reference_source = Source.T(optional=True)
    sources = List.T(Source.T(), optional=True) 
    targets = List.T(Target.T(), optional=True) 
    engine = Engine.T(optional=True)
    store_id = String.T(optional=True)

    depths = List.T()
    misfit_setup = trace.MisfitSetup.T()
    source_time_function = List.T(List.T(Float.T()))
    number_of_time_shifts = Int.T()
    percentage_of_shift = Float.T()
    phase_ids_start = String.T()
    channel_map = Dict.T(String.T(), Int.T(), 
                         optional=True, 
                         default={'N':1, 'E':2, 'Z':3})

    # duration of stencil: 
    marker_perc_length = Float.T(default=1.0)

    # time shift of stencil:
    marker_shift_frac = Float.T(default=0.3)

    test_parameter = String.T(optional=True, default=None)
    test_parameter_value = Float.T(optional=True, default=None)


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
    misfits = Dict.T(Source.T(), Float.T())
