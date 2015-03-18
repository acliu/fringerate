#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C, scipy
from mpl_toolkits.basemap import Basemap

FQ = .159
aa = a.cal.get_aa('psa898_v003', n.array([FQ]))

h = a.healpix.HealpixMap(nside=64)
def beam_area(dat):
    return n.sum(dat) * 4 * n.pi / h.npix()

xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
tx,ty,tz = n.dot(aa._eq2zen, xyz)
_bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]
_bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]
bmXX = n.where(tz > 0, _bmx**2, 0)
bmYY = n.where(tz > 0, _bmy**2, 0)

bl = aa.get_baseline(0,16,'r') * FQ
print 'Baseline:', bl
fng = C.frf_conv.mk_fng(bl, *xyz)
binwidth = .00005
bin_edges = n.arange(-.01+binwidth/2,.01,binwidth)
hXX, bin_edges = n.histogram(fng, bins=bin_edges, weights=bmXX*bmXX)
hXY, bin_edges = n.histogram(fng, bins=bin_edges, weights=bmXX*bmYY)
hYY, bin_edges = n.histogram(fng, bins=bin_edges, weights=bmYY*bmYY)
bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
#polmatch = n.where(hYY > 0, hXY / hYY, 0)
polmatch = n.where(hYY > 0, n.sqrt(hXX / hYY), 0)
polmatch = scipy.interpolate.interp1d(bins, polmatch, kind='linear', bounds_error=False, fill_value=0)
hXX = n.sqrt(hXX)
hXY = n.sqrt(hXY)
hYY = n.sqrt(hYY)

hXX /= hXX.max()
hYY /= hYY.max()
#p.plot(bins, hXX,'k'); p.plot(bins, hYY,'k')#; p.show()
#p.plot(bins, n.where(hXX>0,hYY/hXX,0)); p.show()
bm_pol = polmatch(fng)
bmXXm = bmXX 
bmYYm = bmYY * bm_pol

hpol, bin_edges = n.histogram(fng, bins=bin_edges, weights=(bmXXm-bmYYm)**2)
htot, bin_edges = n.histogram(fng, bins=bin_edges, weights=(bmXXm+bmYYm)**2)
bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
match = n.where(hpol > 0, n.sqrt(htot/hpol), 0); match /= match.max() # XXX sqrt or no?
wgt_pol = scipy.interpolate.interp1d(bins, match, kind='linear', bounds_error=False, fill_value=0)
tbins,firs,frbins,frfs= C.frf_conv.get_fringe_rate_kernels(n.array([wgt_pol]),42.9,403)
frfs[0] /= frfs[0].max()
wgt_pol = scipy.interpolate.interp1d(frbins, frfs[0], kind='linear')


frf,bins,wgt,(cen,wid) = C.frf_conv.get_optimal_kernel_at_ref(aa, 0, (0,16))

bwfrs = C.frf_conv.get_beam_w_fr(aa, (0,16), ref_chan=0)
tbins,firs,frbins,frfs = C.frf_conv.get_fringe_rate_kernels(bwfrs,42.9,403)
wgt = scipy.interpolate.interp1d(frbins, frfs[0], kind='linear')

wgt_bins = wgt(bins); wgt_bins /= wgt_bins.max()
wgt_pol_bins = wgt_pol(bins); wgt_pol_bins /= wgt_pol_bins.max()
wgt_both_bins = wgt_bins * wgt_pol_bins; wgt_both_bins /= wgt_both_bins.max()

fig = p.figure(figsize=(6,6))
fig.subplots_adjust(left=.15, top=.95, bottom=.12, wspace=.15, hspace=.24, right=0.95)

p.subplot(211)
p.plot(1e3*bins,frf,'k')
p.plot(1e3*bins,wgt_bins,'r')
p.plot(1e3*bins,polmatch(bins),'b--')
p.plot(1e3*bins,wgt_pol_bins,'b')
p.plot(1e3*bins,wgt_both_bins,'g')
p.grid()
p.ylim(-0.1,2.)
p.xlim(-1,1.5)
p.xlabel('Fringe Rate [mHz]')
p.ylabel('Filter Response')

firs[0] /= n.abs(firs[0]).max()
p.subplot(212)
p.plot(tbins, n.abs(firs[0]), '#7f3f0e')
p.plot(tbins, n.real(firs[0]), 'c')
p.plot(tbins, n.imag(firs[0]), 'm')
p.ylim(-1.1,1.1)
p.xlim(-5000,5000)
p.xlabel('Time [s]')
p.ylabel('Filter Response')
p.grid()
p.show()

pol_wgt = wgt_pol(fng).real
fng_wgt = wgt(fng).real
#bmXXm_fng = bmXXm * fng_wgt
#bmYYm_fng = bmYYm * fng_wgt

bmI = 0.5 * (bmXX + bmYY)
bmQ = 0.5 * (bmXX - bmYY)
bmIm = 0.5 * (bmXXm + bmYYm)
bmQm = 0.5 * (bmXXm - bmYYm)

bmIm_fng = bmIm * fng_wgt
bmQm_fng = bmQm * fng_wgt
bmIm_fngpol = bmIm_fng * pol_wgt
bmQm_fngpol = bmQm_fng * pol_wgt

print 'Unweighted:', n.sqrt(n.sum(bmQ**2)/n.sum(bmI**2))
print 'Matched:   ', n.sqrt(n.sum(bmQm**2)/n.sum(bmIm**2))
print 'Optimal SN:', n.sqrt(n.sum(bmQm_fng**2)/n.sum(bmIm_fng**2))
print 'Opt pol+SN:', n.sqrt(n.sum(bmQm_fngpol**2)/n.sum(bmIm_fngpol**2))

print 'Unweighted:', (n.sum(bmQ**2)/n.sum(bmI**2))
print 'Matched:   ', (n.sum(bmQm**2)/n.sum(bmIm**2))
print 'Optimal SN:', (n.sum(bmQm_fng**2)/n.sum(bmIm_fng**2))
print 'Opt pol+SN:', (n.sum(bmQm_fngpol**2)/n.sum(bmIm_fngpol**2))

print 'Whole PP beam area:', beam_area(bmI**2)
print 'Matched PP beam area:', beam_area(bmIm**2)
print 'Optimal PP beam area:', beam_area(bmIm_fng**2)
print 'Opt pol PP beam area:', beam_area(bmIm_fngpol**2)

