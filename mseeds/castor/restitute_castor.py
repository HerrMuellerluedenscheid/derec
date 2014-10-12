from pyrocko import pz, io, trace, rdseed, pile
from collections import defaultdict
import glob

pzs = defaultdict()
seeds= defaultdict()

pz_fns = glob.glob('RESPONSES/PAZ/*')
print 'read pz'
for fn in pz_fns:
    fn_splitted = fn.split('_', 5)[2:5]
    fn_splitted.insert(2,'')
    zeros, poles, constant = pz.read_sac_zpk(fn)
    zeros.append(0.0j)
    rest_sts2 = trace.PoleZeroResponse(poles, zeros, 1./constant)
    pzs[tuple(fn_splitted[1:])] = rest_sts2

#fns_traces = glob.glob('2013-10-02T23-06-50/2013*')
fns_traces = glob.glob('2013-10-01T03-32-45_new/2013*')
#fns_traces = glob.glob('2013-10-02T23-29-31/2013*')
p = pile.make_pile(fns_traces)
traces = []
for fn in fns_traces:
    traces.extend(io.load(fn))


rdseed_files = glob.glob('RESPONSES/RESP*')
for rdseed_file in rdseed_files:
    stat = rdseed_file.split('.')
    seedvolume = rdseed.SeedVolumeAccess(rdseed_file, datapile=p)
    seeds[tuple(stat[1:])] = seedvolume


out_traces = []

for tr in traces:

    # create pole-zero response function object for restitution, so poles and zeros
    # from the response file are swapped here.
    try:
        rest_sts2=pzs[tr.nslc_id[1:]]
    except KeyError:
        rest_sts2=seeds[tr.nslc_id].get_restitution(tr,'evalresp')

    tr.ydata -= tr.ydata.mean()
    displacement =  tr.transfer(
        10.,                       # rise and fall of time domain taper in [s]
        (0.05, 0.1, 20., 40.),     # frequency domain taper in [Hz]
        transfer_function=rest_sts2)

    # change channel id, so we can distinguish the traces in a trace viewer.
    displacement.set_codes(channel=tr.channel[-1])

    out_traces.append(displacement)
io.save(out_traces, 'displacement2.mseed')
