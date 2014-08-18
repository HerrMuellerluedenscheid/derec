from collections import defaultdict
from pyrocko import cake
from pyrocko.gui_util import PhaseMarker
from pyrocko import gf


def cake_first_arrival(distance, depth, model, phase_ids=None):
    """
    Get very first arrival of *phases*. 
    """
    phases = []
    for pid in phase_ids:
        try:
            phases.append(cake.PhaseDef(pid))
        except cake.PhaseDefParseError:
            continue

    arrivals = model.arrivals(distances=[distance*cake.m2d], 
                            phases=phases,
                            zstart=depth)

    try:
        tmin = min(arrivals, key=lambda x: x.t).t
    except ValueError:
        print 'cake says: None of defined phases availble in that area. Try another phase'
        raise

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
        self.as_dict = defaultdict(dict)

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
            if not use_cake:
                if isinstance(self.phase_ids_start, list):
                    phase_ids_start_fomosto = '|'.join(self.phase_ids_start)
                else:
                    phase_ids_start_fomosto = self.phase_ids_start

                try:
                    tmin = self.store.t('first(%s)'% phase_ids_start_fomosto, key)
                except gf.store.NoSuchPhase as e:
                    print str(e), ' using cake instead'
                    tmin = None

            if use_cake or tmin==None:
                tmin = cake_first_arrival(dist, source.depth, self.model,
                        phase_ids=self.phase_ids_start)

            tmin += source.time
            self.tmin_phase_cache[key] = tmin

        m = PhaseMarker(nslc_ids=[target.codes],
                        tmin=tmin,
                        tmax=tmin,
                        kind=1,
                        event=None,
                        phasename='cached')

        self.as_dict[source][target] = m

        if self.tmax_phase_cache.get(key, False):
            tmax = self.tmax_phase_cache[key]
        else:
            if perc is not None:
                tmax = tmin + static_length + (tmin-source.time) * perc / 100.

            elif perc is None and use_cake:
                print 'use cake'
                tmax = cake_first_arrival(dist, 
                                          source.depth, 
                                          self.model,
                                          phases=self.phase_ids_end)

            else:
                print 'use fomosto'
                tmax = store.t('first(%s)'%self.phase_ids_end, key)
                tmax += source.time

            self.tmax_phase_cache[key] = tmax

        return tmin, tmax

    def fill_cache(self, sources, targets, **kwargs):
        for s in sources:
            for t in targets:
                self.get_cached_arrivals(target=t,source=s, **kwargs)


    def __str__(self):
        st = 'phase_ids_start: '
        for id in self.phase_ids_start:
            st+= '%s| '%id
        st+='\n'
        
        for k in self.tmin_phase_cache.keys()[0:9]:
            st+= str(k)+str(self.tmin_phase_cache[k])

        return st


