#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys, glob

#files = glob.glob('lst.2455903.[3-5]*.uvAs')
files = glob.glob('lst.2455903.*.uvAs')

dsum,dwgt = {}, {}
for f in files:
    print 'Reading', f, f+'F'
    t,d1,f1 = C.arp.get_dict_of_uv_data([f],'0_16', 'I')
    t,d2,f2 = C.arp.get_dict_of_uv_data([f+'F'],'0_16', 'I')
    for bl in d1:
        dsum[bl], dwgt[bl] = 0., 0.
        for pol in d1[bl]:
            w = n.where(f2[bl][pol], 0, n.where(n.abs(d2[bl][pol]) == 0, 0, 1))
            dsum[bl] += n.sum(n.abs(w*d2[bl][pol])**2, axis=0)
            dwgt[bl] += n.sum(n.abs(w*d1[bl][pol])**2, axis=0)
for bl in dsum:
    p.plot(n.sqrt(dsum[bl]/dwgt[bl]), label=f+str(a.miriad.bl2ij(bl)))
    ans = n.sqrt(dsum[bl].sum()/dwgt[bl].sum())
    print 'Signal Gain (Jy):', ans
    print 'Relative Beam Area:', ans**2
p.show()
        
            
