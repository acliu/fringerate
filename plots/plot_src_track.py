#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import glob

def get_data(globstr, ch=100):
    d = glob.glob(globstr)
    t, d, _ = C.arp.get_dict_of_uv_data(d, 'all', 'I')
    t = n.arange(t.size) * (t[1]-t[0]) / a.ephem.second
    d = d[273]['I'][:,ch]
    return t,d

ax1 = p.subplot(311)
t00, d00 = get_data('data/flat_beam/d+00/*uvAs')
p.plot(t00, d00.real, 'k')
t00, d00 = get_data('data/flat_beam/d+00/*uvAsF')
p.plot(t00, d00.real, 'r')
p.ylim(-1.5,1.5)
p.setp(ax1.get_xticklabels(), visible=False)

ax2 = p.subplot(312, sharex=ax1)
t30, d30 = get_data('data/flat_beam/d-30/*uvAs')
p.plot(t30, d30.real, 'k')
t30, d30 = get_data('data/flat_beam/d-30/*uvAsF')
p.plot(t30, d30.real, 'r')
p.ylim(-1.5,1.5)
p.setp(ax2.get_xticklabels(), visible=False)

ax3 = p.subplot(313, sharex=ax1)
t60, d60 = get_data('data/flat_beam/d-60/*uvAs')
p.plot(t60, d60.real, 'k')
t60, d60 = get_data('data/flat_beam/d-60/*uvAsF')
p.plot(t60, d60.real, 'r')
p.ylim(-1.5,1.5)

p.show()
