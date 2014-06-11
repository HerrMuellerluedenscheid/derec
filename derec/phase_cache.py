from collections import defaultdict
from pyrocko import cake


def cake_first_arrival(distance, depth, model, phase_ids=None):
    """
    Get very first arrival of *phases*. 
    """
    phases = [cake.PhaseDef(ph) for ph in phase_ids]
    arrivals = model.arrivals([distance*cake.m2d], 
                            phases,
                            zstart=depth,
                            zstop=depth)

    tmin = min(arrivals, key=lambda x: x.t).t
    return tmin



class PhaseCache():
    def __init__(self, tmin_phase_cache=None, tmax_phase_cache=None, store=None,
            phase_ids_start=None, phase_ids_end=None):

        self.tmin_phase_cache = tmin_phase_cache if tmin_phase_cache else defaultdict()
        self.tmax_phase_cache = tmax_phase_cache if tmax_phase_cache else defaultdict()
        self.model = store.config.earthmodel_1d
        self.store = store
        self.phase_ids_start = phase_ids_start if phase_ids_start else ['p','P']
        self.phase_ids_end = phase_ids_end 

    def flush(self):
        self.tmin_phase_cache = defaultdict()
        self.tmax_phase_cache = defaultdict()

    def get_cached_tmin(self, target, source):
        dist = source.distance_to(target)
        key = (source.depth, dist)
        return self.tmin_phase_cache[key]
      
    def get_cached_arrivals(self, target, source, static_length=0., 
            perc=None, use_cake=True):

        dist = source.distance_to(target)
        key = (source.depth, dist)
        
        try:
            tmin = self.get_cached_tmin(target, source)
        except KeyError:
            if use_cake:
                tmin = cake_first_arrival(dist, source.depth, self.model,
                        phase_ids=self.phase_ids_start)
            else:
                tmin = self.store.t('first(%s)'% self.phase_ids_start, key)
            self.tmin_phase_cache[key] = tmin

        if self.tmax_phase_cache.get(key, False):
            tmax = self.tmax_phase_cache[key]
        else:
            if perc:
                tmax = tmin + static_length + tmin * perc / 100.

            elif use_cake:
                print 'use cake'
                tmax = cake_first_arrival(dist, source.depth, self.model,
                        phases=self.phase_ids_end.split('|'))

            else:
                print 'use fomosto'
                tmax = store.t('first(%s)'%self.phase_ids_end, key)
            self.tmax_phase_cache[key] = tmax


        tmin += source.time
        tmax += source.time

        return tmin, tmax


