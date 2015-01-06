#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import glob

def get_data(globstr, antstr='all', polstr='I', decimate=1, decphs=0, verbose=False, recast_as_array=True, ch=100):
    filenames = glob.glob(globstr)
    times, dat = [], {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        if decimate > 1: uv.select('decimate', decimate, decphs)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            t = uv['lst']
            if len(times) == 0 or t != times[-1]: times.append(t)
            bl = a.miriad.ij2bl(i,j)
            if not dat.has_key(bl): dat[bl] = {}
            pol = a.miriad.pol2str[uv['pol']]
            if not dat[bl].has_key(pol):
                dat[bl][pol] = []
            dat[bl][pol].append(d)
    if recast_as_array:
        for bl in dat.keys():
          for pol in dat[bl].keys():
            dat[bl][pol] = n.array(dat[bl][pol])
    times = n.array(times)
    times = n.where(times > 5, times-2*n.pi,times)
    return times-n.pi/2, dat.values()[0][polstr][:,ch]

fig = p.figure(figsize=(5,5))
fig.subplots_adjust(left=.15, top=.95, bottom=.15, wspace=.15, hspace=.1, right=0.95)
ax1 = p.subplot(311)
t00, d00 = get_data('data/flat_beam/d+00/*uvAs')
print t00.shape
p.plot(t00, d00.real, 'k')
t00, d00 = get_data('data/flat_beam/d+00/*uvAsF')
p.plot(t00, d00.real, 'r')
p.ylim(-1.5,1.5)
p.yticks([-1,0,1])
p.xlim(-1.3,1.3)
p.setp(ax1.get_xticklabels(), visible=False)

ax2 = p.subplot(312, sharex=ax1)
t30, d30 = get_data('data/flat_beam/d-30/*uvAs')
p.plot(t30, d30.real, 'k')
t30, d30 = get_data('data/flat_beam/d-30/*uvAsF')
p.plot(t30, d30.real, 'r')
p.ylim(-1.5,1.5)
p.yticks([-1,0,1])
p.xlim(-1.3,1.3)
p.ylabel('Normalized Fringe Amplitude')
p.setp(ax2.get_xticklabels(), visible=False)

ax3 = p.subplot(313, sharex=ax1)
t60, d60 = get_data('data/flat_beam/d-60/*uvAs')
p.plot(t60, d60.real, 'k')
t60, d60 = get_data('data/flat_beam/d-60/*uvAsF')
p.plot(t60, d60.real, 'r')
p.ylim(-1.5,1.5)
p.yticks([-1,0,1])
p.xlim(-1.3,1.3)
p.xlabel('Hour Angle (radians)')

p.show()
