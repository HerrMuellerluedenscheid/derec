from pyrocko import trace, io
import matplotlib.pyplot as plt
import numpy as num
import os
pjoin = os.path.join
'''
restitude seismic traces using resp files. 

Give the number of stations, that you would like to restitude as
*stations_needed*.
'''
def save_noise(directory, traces):
    if not os.path.exists(directory):
        os.makedirs(directory)

    for t in traces:
        io.save(traces, pjoin(directory,'%s.%s.%s.%s'%(t.nslc_id)))


def equalize_number_of_samples(traces):
    minlen = min([len(t.get_ydata()) for t in traces])
    map(lambda t: t.set_ydata(t.get_ydata()[:minlen]), traces)
        

def check_noise(traces):
    verticals = filter(lambda x: x.nslc_id[3]=='BHZ', traces)
    horizontals = list(set(traces)-set(verticals))

    checked = check_ampspecs(verticals)
    checked.extend(check_ampspecs(horizontals))

    plt.show()
    return checked

def check_ampspecs(traces, **kwargs):
    ampspecs = []
    checked = []
    log_ydata = {}
    for t in traces:
        fx,fa = t.spectrum(tfade=120)
        abs_fa = num.abs(fa)
        log_abs_fa = num.log(abs_fa)
        log_ydata.update({t:log_abs_fa})
    log_ydata_array = num.array(log_ydata.values())
    log_median = num.median(log_ydata_array, axis=0)
    sigma = num.std(log_median)
    mean_log_median = num.mean(log_median)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for t, logdata in log_ydata.items():
        mean_log_y = num.mean(logdata)
        if num.abs(mean_log_median-mean_log_y)<0.3333*sigma:
            alpha = 1.
            symbol = 'x'
            checked.append(t)
        else:
            alpha = 0.05
            symbol = '.'
        
        ax.plot(fx, num.exp(logdata), symbol, alpha=alpha)
                
    ax.plot(fx, num.exp(log_median),'--')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.draw()
    return checked
 


traces = io.load('/scratch/local1/marius/GSN.mseed')

channels = ['BHE', 'BHN', 'BHZ']
respdir = 'GSN_resp/'
restdir = 'restitute/'
stations_needed = 40
i=0

done = []
c_left = set(channels)
#t_stack = []
t_stack = {}
ready = []
for t_i, t in enumerate(traces):
    if not t.channel in channels:
        continue
    else:
        try:
            t_stack[t.nslc_id[:3]].append(t)
            if len(t_stack[t.nslc_id[:3]])==3:
                print '-----'
                for t in t_stack[t.nslc_id[:3]]:
                    t.downsample_to(1.0)
                    tcodes = t.nslc_id
                    invevalresp = trace.InverseEvalresp(respdir, t)
                    print 'restitute: %s'%t
                    tnew = t.transfer(tfade=120., 
                                      freqlimits=[0.01, 0.02, 10, 20], 
                                      transfer_function=invevalresp)

                    done.append(tnew)
                i+=1
                if i==stations_needed:
                    break
        except KeyError:
            t_stack.update({t.nslc_id[:3]: [t]})

equalize_number_of_samples(done)
checked = check_noise(done)

save_noise('checked_noise', checked)



