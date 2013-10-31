from pyrocko import trace, util, io
import numpy as num

nsamples = 100
tmin = util.str_to_time('2010-02-20 15:15:30.100')
data1 = num.random.random(nsamples)
data2 = num.random.random(nsamples)
t1 = trace.Trace(network='test', station='TEST', channel='Z', deltat=0.5, tmin=tmin, ydata=data1)
t2 = trace.Trace(network='test', station='TEST', channel='N', deltat=0.5, tmin=tmin, ydata=data2)
io.save([t1,t2], 'my_precious_test_traces.mseed')            # all traces in one file

