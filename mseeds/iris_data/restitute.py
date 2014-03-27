from pyrocko import trace, io, InverseEvalresp

traces = io.load('GSN.mseed')

channels = ['BHE', 'BHN', 'BHZ']
respdir = 'GSN_resp'
stations_needed = 13
for i in stations_needed:
    for t in traces:
        tcodes = t.nslc_codes
        respfile = respdir + 'RESP.' + tcodes 
        invevalresp = trace.InverseEvalresp(respfile)
        t.transfer(tfade=10., freqlimits=[], transfer_function=invevalresp)
            
        

