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
i=1

done = [[]]
for t in traces:
    if not all(map(lambda x: done.count(x)==stations_needed, channels)):    
        if t.channel in channels:
            if done.count(t.nslc_id[3])== stations_needed:
                continue
            tcodes = t.nslc_id
            respfile = respdir + 'RESP.%s.%s.%s.%s.'%tcodes 
            invevalresp = trace.InverseEvalresp(respdir, t)
            print 'restitute: %s'%t
            tnew = t.transfer(tfade=10., 
                              freqlimits=[0.01, 0.02, 10, 20], 
                              transfer_function=invevalresp)

            io.save(tnew, restdir+'%s.%s.%s.%s'%(t.nslc_id))
            done.append(tnew.nslc_id[3])
            i+=1
            
