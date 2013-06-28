#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import scipy.interpolate

FQ = .164
#BINWIDTH = 1. / 7200. # Hz
#BINWIDTH = 1. / 3600. # Hz
BINWIDTH = .00005
bin_edges = n.arange(-.01+BINWIDTH/2,.01,BINWIDTH)
NBINS = len(bin_edges)

def mk_fng(bl, ex, ey, ez):
    #f = -(bl[0]*n.sin(lons) + bl[1]*n.cos(lons)) * n.cos(lats)
    return 2*n.pi/a.const.sidereal_day * (bl[0]*ex + bl[1]*ey) * n.sqrt(1 - ez**2) # Hz

def bm_area(hmap): return n.sum(hmap) * 4*n.pi / hmap.size
def bm2_area(hmap): return n.sum(hmap**2) * 4*n.pi / hmap.size

aa = a.cal.get_aa('psa898_v003', n.array([FQ]))
if False:
    for ant in aa: ant.set_pointing(twist=n.pi/4)
if False: aa.lat = '0:00'; aa.update()
lat = aa.lat
print lat
h = a.healpix.HealpixMap(nside=64)
DIM = 400
RES = .5
SIZE = DIM / RES
im = a.img.Img(DIM, res=RES)

def proj_hmap(dat, mode='e'):
    h.map = dat
    #if mode == 'e': x,y,z = im.get_eq(dec=aa.lat, center=(SIZE/2,SIZE/2))
    if mode == 'e': x,y,z = im.get_eq(dec=lat, center=(SIZE/2,SIZE/2))
    elif mode == 't': x,y,z = im.get_top(center=(SIZE/2,SIZE/2))
    else: raise ValueError(mode)
    v = n.logical_not(x.mask)
    x,y,z = x.flatten(), y.flatten(), z.flatten()
    d = h[x,y,z]
    d.shape = (SIZE,SIZE)
    d = n.where(v, d, 0)
    return d

def beam_area(dat):
    return n.sum(dat) * 4 * n.pi / h.npix()

xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
tx,ty,tz = n.dot(aa._eq2zen, xyz)
_bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]
_bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]
bmxx = n.where(tz > 0, _bmx**2, 0)
bmyy = n.where(tz > 0, _bmy**2, 0)
bm_I = 0.5 * (bmxx + bmyy)
bm_Q = 0.5 * (bmxx - bmyy)

xyz = (xyz[1],xyz[0],xyz[2])
#bl = n.array([100, 0, 0])
bl = aa.get_baseline(0,16,'r') * FQ
#bl = aa.get_baseline(0,28,'r') * FQ
#bl = aa.get_baseline(0,3,'r') * 10 * FQ
print 'Baseline:', bl
fng = mk_fng(bl, *xyz)
print fng.max(), fng.min()

h_I, bin_edges = n.histogram(fng, bins=bin_edges, weights=bm_I**2) # squared to correspond to EoR sensitivity
hxx, bin_edges = n.histogram(fng, bins=bin_edges, weights=bmxx**2) # squared to correspond to EoR sensitivity
hyy, bin_edges = n.histogram(fng, bins=bin_edges, weights=bmyy**2) # squared to correspond to EoR sensitivity
h_Q, bin_edges = n.histogram(fng, bins=bin_edges, weights=bm_Q**2) # squared to correspond to EoR sensitivity
bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
h_I /= h_I.max()
wgt = scipy.interpolate.interp1d(bins, h_I, kind='linear')
ratio = n.where(hyy == 0, 0, n.sqrt(hxx/hyy))
wgt_ratio = scipy.interpolate.interp1d(bins, ratio, kind='linear')

wgt_ratio = wgt_ratio(fng)
bm_I_ = 0.5 * (bmxx + bmyy * wgt_ratio)
bm_Q_ = 0.5 * (bmxx - bmyy * wgt_ratio)
h_I_, bin_edges = n.histogram(fng, bins=bin_edges, weights=bm_I_**2) # squared to correspond to EoR sensitivity
h_Q_, bin_edges = n.histogram(fng, bins=bin_edges, weights=bm_Q_**2) # squared to correspond to EoR sensitivity
leak = h_I_/h_Q_.clip(1e-6,n.Inf); leak /= leak.max()
wgt_leak = scipy.interpolate.interp1d(bins, leak, kind='linear')

if True:
    #p.subplot(211)
    p.plot(bins, h_I, 'k')
    p.plot(bins, leak, 'r')
    #p.plot(bins, h_I_ * 4*n.pi/h.npix(), 'g')

    #p.subplot(212)
    p.plot(bins, ratio)

    def gauss(bins, cen, wid):
        return n.exp(-(bins-cen)**2/(2*wid**2))
        
    def fit_gaussian(h, prms):
        cen,wid = prms
        g = gauss(bins, cen, wid)
        g = n.where(bins > fng.max(), 0, g)
        score = n.sum((h-g)**2)
        return score

    (cen,wid), score = a.optimize.fmin(lambda prms: fit_gaussian(h_I,prms), [.001, .0001], 
            full_output=1, disp=0, maxfun=1000, maxiter=n.Inf, ftol=1e-6, xtol=1e-6)[:2]
    print 'Stokes I:', cen, wid

    #p.subplot(211)
    p.plot(bins, n.exp(-(bins-cen)**2/(2*wid**2)), 'k:')
    #for sub in [211,212]:
    if True:
        #p.subplot(sub)
        p.xlabel('Fringe Rate (Hz)')
        p.ylabel('Normalized Beam Response')
        p.xlim(-.0005, .0015)
    p.show()

def plot_hmap(dat, mode='log', mx=0, drng=3):
    print bm_area(dat), bm2_area(dat)
    d = proj_hmap(dat,'e')
    if mode == 'log': d = n.log10(n.abs(d).clip(10**(mx-drng),n.Inf))
    p.imshow(d, vmax=mx, vmin=mx-drng, origin='lower')

wgt_snr = wgt(fng)
wgt_leak = wgt_leak(fng)

p.subplot(241)
plot_hmap(bm_I_)

p.subplot(245)
plot_hmap(bm_I_ * wgt_snr)

p.subplot(242)
plot_hmap(bm_Q_)

p.subplot(246)
plot_hmap(bm_Q_ * wgt_snr)

p.subplot(243)
plot_hmap(bm_I_ * wgt_leak)

p.subplot(247)
plot_hmap(bm_I_ * wgt_snr * wgt_leak)

p.subplot(244)
plot_hmap(bm_Q_ * wgt_leak)

p.subplot(248)
plot_hmap(bm_Q_ * wgt_snr * wgt_leak)

#bm_diff = n.abs(bmxx - bmyy)**2
#bm_sum = n.abs(bmxx + bmyy)**2
#leakage = n.sqrt(n.sum(bm_diff) / n.sum(bm_sum))
#print 'Raw leakage:', leakage
#print 'Matched leakage:', n.sqrt(n.sum(n.abs(bmxx - bmyy_)**2)/n.sum(n.abs(bmxx + bmyy_)**2))
#print 'Weighted leakage:', n.sqrt(n.sum(n.abs(wgt_leak*(bmxx - bmyy_))**2)/n.sum(n.abs(wgt_leak*(bmxx + bmyy_))**2))
#
#bm_diff_fng = n.abs(wgt_snr*(bmxx-bmyy_))**2
#leakage_fng = n.sqrt(n.sum(bm_diff_fng)/n.sum(bm_sum))
#print 'Fringe leakage relative to orginal beam area:', leakage_fng
#print 'Fractional leakage, Fng beam:', n.sqrt(n.sum(bm_diff_fng)/n.sum(n.abs(wgt_snr*(bmxx+bmyy_))**2))

p.show()

