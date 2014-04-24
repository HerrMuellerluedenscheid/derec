from pyrocko import trace, io
'''
restitude seismic traces using resp files. 

Give the number of stations, that you would like to restitude as
*stations_needed*.
'''

traces = io.load('GSN.mseed')

channels = ['BHE', 'BHN', 'BHZ']
respdir = 'GSN_resp/'
restdir = 'restitute/'
stations_needed = 13
i=0

done = [[]]
c_left = set(channels)
t_stack = []
for t_i, t in enumerate(traces):
    if not t.channel in channels:
        t_stack = []
        continue
    t_stack.append(t)
    if t.nslc_id[:3]!=t_stack[0].nslc_id[:3]:
        t_stack = [t]
        continue
    if len(t_stack)==3:
        print '-----'
        for t in t_stack:
            tcodes = t.nslc_id
            invevalresp = trace.InverseEvalresp(respdir, t)
            print 'restitute: %s'%t
            tnew = t.transfer(tfade=10., 
                              freqlimits=[0.01, 0.02, 30, 60], 
                              transfer_function=invevalresp)

            io.save(tnew, restdir+'%s.%s.%s.%s'%(t.nslc_id))
            done.append(tnew.nslc_id[3])
        i+=1
        if i==stations_needed:
            break
