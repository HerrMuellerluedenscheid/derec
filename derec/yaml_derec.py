import numpy as num

from pyrocko.guts import Object, Int, Float, List, String, Tuple, Timestamp, Dict
from pyrocko.guts_array import Array
from pyrocko.gf import *
from pyrocko import trace

from derec import derec_utils as du

guts_prefix = 'derec.yaml_derec'

class yamlMarker(Object):
    """
    yaml version of the pyrocko.gui_util.Marker class
    """
    nslc_ids = List.T(Tuple.T(4, String.T()))
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


class TestCaseSetup(Object):
    tagname = String.T(default='TestCaseSetup')
    reference_source = Source.T(optional=True)
    sources = List.T(Source.T(), optional=True, default=None) 
    targets = List.T(Target.T(), optional=True) 
    engine = Engine.T(optional=True)
    store_id = String.T(optional=True)

    # depths have to given as a list of floats!
    depths = List.T(Float.T())
    misfit_setup = trace.MisfitSetup.T()
    source_time_function = List.T(List.T(Float.T()))
    number_of_time_shifts = Int.T(optional=True, default=0)
    percentage_of_shift = Float.T(optional=True)
    time_shift = Float.T(optional=True)
    phase_ids_start = List.T(String.T(optional=True, default=['p','P']))
    channel_map = Dict.T(String.T(), Int.T(), 
                         optional=True, 
                         default={'N':1, 'E':2, 'Z':3})

    static_length = Float.T(optional=True, default=0.00001)

    # to be added to duration of stencil: 
    marker_perc_length = Float.T(default=1.0)

    # time shift of stencil:
    marker_shift_frac = Float.T(default=0.3)

    test_parameter = String.T(optional=True, default=None)
    test_parameter_value = Float.T(optional=True, default=None)

    outlier_threshold = Float.T(optional=True, default=10.)

    #def __init__(self, **kwargs):
    #    self.__dict__.update(kwargs)

    #    for k,v in kwargs.items():
    #        setattr(self, k, v)
    #        print getattr(self, k, v)

    #    @property
    #    def sources(self):
    #        if not self._sources:
    #            print self._sources
    #            return self._sources
    #        else:
    #            self._sources = du.test_event_generator(self.reference_source,\
    #                                               self.depths)
    #            
    #            self._sources.regularize()

    #            return self._sources


class TestCaseData(Object):
    """
    A TestCaseData object contains all waveforms, misfit results as dictionary
    and a TestCaseSetup object.
    """
    description = String.T(optional=True)
    references = Dict.T(Source.T(), Dict.T(Target.T(), SeismosizerTrace.T()), 
            optional=True)
    candidates = Dict.T(Source.T(), Dict.T(Target.T(), SeismosizerTrace.T()), 
            optional=True)
    processed_references = Dict.T(Source.T(), Dict.T(Target.T(), SeismosizerTrace.T()),
            optional=True)
    processed_candidates = Dict.T(Source.T(), Dict.T(Target.T(), SeismosizerTrace.T()), 
            optional=True)

    candidates_markers = Dict.T(Source.T(), Dict.T(Target.T(), yamlMarker.T()))
    reference_markers = Dict.T(Source.T(), Dict.T(Target.T(), yamlMarker.T()))
    test_case_setup = TestCaseSetup.T(optional=True)
    misfits = Dict.T(Source.T(), Float.T())
