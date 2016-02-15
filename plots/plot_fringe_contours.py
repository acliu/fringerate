#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
from mpl_toolkits.basemap import Basemap

FQ = .164
BINWIDTH = 1. / 7200. # Hz
#BINWIDTH = 1. / 3600. # Hz
bin_edges = n.arange(-.01+BINWIDTH/2,.01,BINWIDTH)
NBINS = len(bin_edges)

def mk_fng(bl, ex, ey, ez):
    #f = -(bl[0]*n.sin(lons) + bl[1]*n.cos(lons)) * n.cos(lats)
    return -2*n.pi/a.const.sidereal_day * (bl[0]*ex + bl[1]*ey) * n.sqrt(1 - ez**2) # Hz
from frf_conv import mk_fng

aa = a.cal.get_aa('psa898_v003', n.array([FQ]))
if False:
    for ant in aa: ant.set_pointing(twist=n.pi/4)
if False: aa.lat = '0:00'; aa.update()
lat = aa.lat
print lat
fig = p.figure(figsize=(4,5))
fig.subplots_adjust(left=.05, top=.95, bottom=.05, right=0.95)
m = Basemap(projection='ortho', lat_0=lat, lon_0=180, rsphere=1.)
h = a.healpix.HealpixMap(nside=128)
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
    d = n.where(v, d, n.NaN)
    return d

def beam_area(dat):
    return n.sum(dat) * 4 * n.pi / h.npix()

xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
tx,ty,tz = n.dot(aa._eq2zen, xyz)
if False:
    aa[0].set_pointing(twist=0)
    aa[28].set_pointing(twist=n.pi/4)
_bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]
_bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]
_bma = aa[28].bm_response((tx,ty,tz),pol='x')[0]
_bmb = aa[28].bm_response((tx,ty,tz),pol='y')[0]
#bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]**2
#bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]**2
#bm = bmx - bmy
#bm = _bmx * _bmb #+ _bmy * _bma
#bm = _bmy * _bma
bm = 0.5 * (_bmx**2 + _bmy**2)
bm = n.where(tz > 0, bm, 0)

def nos(dat, t):
    pwr_bm = beam_area(bm)
    pwr2_bm = beam_area(dat**2)
    print pwr2_bm, pwr_bm, t
    return pwr_bm**2 / pwr2_bm / t

#print beam_area(bm), beam_area(bm**2), nos(bm, 1.), nos(bm,1.) / n.sqrt(NBINS)

#xyz = (xyz[1],xyz[0],xyz[2])
#bl = n.array([100, 0, 0])
bl = aa.get_baseline(0,16,'r') * FQ
#bl = aa.get_baseline(0,28,'r') * FQ
#bl = aa.get_baseline(0,3,'r') * 10 * FQ
print 'Baseline:', bl
fng = mk_fng(bl, *xyz)
#print fng.max(), fng.min()

hist, bin_edges = n.histogram(fng, bins=bin_edges, weights=bm**2)
#hist, bin_edges = n.histogram(fng, range=(-.01,.01), bins=NBINS, weights=bm**2) # or bm**2?
f = 0
for cnt,b1 in enumerate(bin_edges[:-1]):
    b2 = bin_edges[cnt+1]
    if cnt % 2 == 0: f = n.where(n.logical_and(fng >= b1, fng < b2), 1., f)

#p.subplot(231)
#p.imshow(proj_hmap(bm,'e'), origin='lower')

#p.subplot(232)
#p.imshow(proj_hmap(fng,'e'), vmax=0, vmin=-.007, origin='lower')
m.imshow(1e3*proj_hmap(fng,'e'), origin='lower')
cb = p.colorbar(shrink=1., orientation='horizontal', pad=.05)
cb.set_label('Fringe Rate [mHz]')

#p.subplot(233)
#p.imshow(proj_hmap(f*bm,'e'), origin='lower', alpha=.5)
m.imshow(proj_hmap(f,'e'), origin='lower', cmap='gist_yarg', alpha=.1)

#p.subplot(234)
##hist, bins = n.histogram(fng, range=(-.01,.01), bins=NBINS, weights=bm) # or bm**2?
#bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
#hist1 = hist/n.sum(hist)
#print 'Hist:', hist1.max()
#p.plot(bins, hist1/hist1.max())

noise_sum, noise_wgt = 0, 0
wgts = []
f = 0
for cnt,b1 in enumerate(bin_edges[:-1]):
    b2 = bin_edges[cnt+1]
    fng_wgt = n.where(n.logical_and(fng >= b1, fng < b2), 1., 0)
    fng_amp = bm * fng_wgt
    _n = nos(fng_amp,NBINS)
    w = 1. / _n**2
    wgts.append(w)
    if w > 0:
        #print cnt, b1, _n, w
        noise_sum += (_n*w)**2
        noise_wgt += w
    f = n.where(n.logical_and(fng >= b1, fng < b2), w, f)
wgts = n.array(wgts)

#p.plot(bins, wgts/wgts.max())

fng_amp = bm * f
print fng_amp.max()
fng_amp /= fng_amp.max()
print 'Whole beam area:', beam_area(bm), beam_area(bm**2)
print 'Filtered beam area:', beam_area(fng_amp), beam_area(fng_amp**2)
print n.sqrt(noise_sum) / noise_wgt

print 'Improvement:', (n.sqrt(noise_sum)/noise_wgt) / (nos(bm,1.) / n.sqrt(NBINS))
    
#p.subplot(235)
#p.imshow(proj_hmap(f,'e'), origin='lower')
#
#p.subplot(236)
#p.imshow(proj_hmap(fng_amp,'e'), origin='lower')
m.drawmapboundary(linewidth=2)
m.drawmeridians(n.arange(0, 360, 30), linewidth=2)
m.drawparallels(n.arange(-90,90,30), linewidth=2)

p.show()
