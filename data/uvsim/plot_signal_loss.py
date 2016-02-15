#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys, glob, scipy



#files = glob.glob('lst.2455903.[3-5]*.uvAs')
files = glob.glob('lst.2455903.[3-6]*.uvAs')
#files = glob.glob('lst.2455903.*.uvAs')

dsum,dwgt = {}, {}
for f in files:
    print 'Reading', f, f+'F'
    t,d1,f1 = C.arp.get_dict_of_uv_data([f],'0_16', 'I')
    #t,d2,f2 = C.arp.get_dict_of_uv_data([f+'F'],'0_16', 'I')
    t,d2,f2 = C.arp.get_dict_of_uv_data([f+'L'],'0_16', 'I')
    for bl in d1:
        dsum[bl], dwgt[bl] = 0., 0.
        for pol in d1[bl]:
            w = n.where(f2[bl][pol], 0, n.where(n.abs(d2[bl][pol]) == 0, 0, 1))
            dsum[bl] += n.sum(n.abs(w*d2[bl][pol])**2, axis=0)
            dwgt[bl] += n.sum(n.abs(w*d1[bl][pol])**2, axis=0)

FQ = .151
DIM = 400
RES = .5
SIZE = DIM / RES

aa = a.cal.get_aa('psa898_v003', n.array([FQ]))
bl = aa.get_baseline(0,16,'r') * FQ
print 'Baseline:', bl

h = a.healpix.HealpixMap(nside=64)
def beam_area(dat): return n.sum(dat) * 4 * n.pi / h.npix()
xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
tx,ty,tz = n.dot(aa._eq2zen, xyz)
_bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]
_bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]
bmI = 0.5 * (_bmx**2 + _bmy**2)
bmI = n.where(tz > 0, bmI, 0)

fng = C.frf_conv.mk_fng(bl, *xyz)
frf,bins,wgt,(cen,wid) = C.frf_conv.get_optimal_kernel_at_ref(aa, 0, (0,16))
bwfrs = C.frf_conv.get_beam_w_fr(aa, (0,16), ref_chan=0)
tbins,firs,frbins,frfs = C.frf_conv.get_fringe_rate_kernels(bwfrs,42.9,403)

frfs[0] = n.fft.fftshift(n.fft.fft(n.fft.ifftshift(firs[0]), axis=-1))
wgt = scipy.interpolate.interp1d(frbins, frfs[0].real, kind='linear')
print 'Beam Area (Power):', beam_area(bmI)
print 'Beam Area (Power^2):', beam_area(bmI**2)
fng_wgt = wgt(fng)
fng_bm = fng_wgt * bmI
skypass = n.where(n.abs(frbins) < .0012, 1., 0)
crit_noise_lev = n.sum(skypass**2) / frbins.size
fng_noise_lev = n.sum(n.abs(frfs[0])**2)/ frbins.size
#p.plot(frbins, skypass)
#p.plot(frbins, n.abs(frfs[0]))
#p.show()
print 'Eff Int Time:', 42.9/crit_noise_lev, 42.9/fng_noise_lev
print 'Beam Area Fringe (Power):', beam_area(fng_bm)
print 'Beam Area Fringe (Power^2):', beam_area(fng_bm**2)
print 'Signal Loss:', n.sum(fng_bm**2)/n.sum(bmI**2)
print 'Noise Redux:', fng_noise_lev / crit_noise_lev
print 'Sensitivity:', n.sqrt(crit_noise_lev / fng_noise_lev) * n.sum(fng_bm**2)/n.sum(bmI**2)

for bl in dsum:
    #p.plot(n.sqrt(dsum[bl]/dwgt[bl]), label=f+str(a.miriad.bl2ij(bl)))
    p.plot(dsum[bl]/dwgt[bl], label=f+str(a.miriad.bl2ij(bl)))
    print n.std(dsum[bl]/dwgt[bl])
    ans = n.sqrt(dsum[bl].sum()/dwgt[bl].sum())
    print 'Signal Gain (Jy):', ans
    print 'Relative Beam Area:', ans**2
p.show()
        
            
