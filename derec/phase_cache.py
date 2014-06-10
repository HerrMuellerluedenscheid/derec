from collections import defaultdict
from pyrocko import cake


def cake_first_arrival(distance, depth, model, phase_ids=None):
    """
    Get very first arrival of *phases*. 
    """
    if not phase_ids:
        phase_ids= ['p','P', 'pP']
    
    phases = [cake.PhaseDef(ph) for ph in phase_ids]
    arrivals = model.arrivals([distance*cake.m2d], 
                            phases,
                            zstart=depth,
                            zstop=depth)

    tmin = min(arrivals, key=lambda x: x.t).t
    return tmin



class PhaseCache():
    def __init__(self, tmin_phase_cache={}, store=None, phase_ids_start=[], 
            phase_ids_end=[]):

        self.tmin_phase_cache = tmin_phase_cache
        #print 'new, tmin_phase_cache', self.tmin_phase_cache
        self.model = store.config.earthmodel_1d
        self.store = store
        self.phase_ids_start = phase_ids_start
        self.phase_ids_end = phase_ids_end

    def flush(self):
        self.tmin_phase_cache = defaultdict()
        #self.model = store.config.earthmodel_1d
        #self.store = store
        #self.phase_ids_start = phase_ids_start
        #self.phase_ids_end = phase_ids_end
        print 'phase cache: "flushed myself..." ' , self.tmin_phase_cache

    def get_cached_tmin(self, target, source):
        dist = source.distance_to(target)
        key = (source.depth, dist)
        return self.tmin_phase_cache[key]
      
    def get_cached_arrivals(self, target, source, static_length=0., 
            perc=0, use_cake=True):
        print static_length, 'gca'

        dist = source.distance_to(target)
        key = (source.depth, dist)
        
        try:
            tmin = self.get_cached_tmin(target, source)
            print 'cool, I got a cached tmin', tmin
        except KeyError:
            if use_cake:
                tmin = cake_first_arrival(dist, source.depth, self.model,
                        phase_ids=self.phase_ids_start)
                print 'mmm, I didnt get tmin-> cak-> cake',  tmin
            else:
                tmin = self.store.t('first(%s)'% self.phase_ids_start, key)
            self.tmin_phase_cache[key] = tmin

        if perc or static_length:
            print static_length,
            tmax = tmin + static_length + tmin * perc / 100.
            print tmax-tmin, 'tmax-tmin'

        elif use_cake:
            print 'use cake'
            tmax = cake_first_arrival(dist, source.depth, self.model,
                    phases=self.phase_ids_end.split('|'))

        else:
            print 'use fomosto'
            tmax = store.t('first(%s)'%self.phase_ids_end, key)


        tmin += source.time
        tmax += source.time

        return tmin, tmax


