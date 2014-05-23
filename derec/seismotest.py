from pyrocko import gf
from pyrocko import trace
from pyrocko import orthodrome
import sys
import os

pjoin = os.path.join

tdict = {'depth':[5e3, 10e3, 15e3], 
         'lat':[11.0,11.1,11.2]}

#sources = [
#   gf.MTSource(
#   lat=10.0,
#   lon=lat,
#   depth=depth,
#   mnn=1.,
#   mee=0.3,
#   mdd=0.1, 
#   mne=0.,
#   med=0.,
#   mnd=0.) for depth in [5e3, 10e3, 15e3] for lat in [11.0,11.1,11.2]]
print orthodrome.distance_accurate50m_numpy(10., 10., 10.5,10.)

sources = [
   gf.DCSource(
   lat=10.0,
   lon=10.0,
   depth=depth,
   magnitude=magnitude,
   strike=33.,
   dip=80.,
   rake=5.) for depth in [5e3] for magnitude in [2.]]

targets = [
    gf.Target(
    codes=('', 'STA', '', component),
    lat=10.5,
    lon=10.0,
    store_id='local1',
    ) for component in 'NEZ']


direc = pjoin(os.environ['DEREC_HOME'], 'fomostos')
engine = gf.LocalEngine(store_superdirs=[direc])
response = engine.process(
            sources=sources,
            targets=targets)

#print response

trace.snuffle(response.pyrocko_traces())
