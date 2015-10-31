#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
from mpl_toolkits.basemap import Basemap

FQ = .164
aa = a.cal.get_aa('psa898_v003', n.array([FQ]))
lat = aa.lat
h = a.healpix.HealpixMap(nside=512)
DIM = 800
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

xyz = h.px2crd(n.arange(h.map.size), ncrd=3)
ty,tx,tz = n.dot(aa._eq2zen, xyz)
xyz = (xyz[1],xyz[0],xyz[2])
_bmx = aa[0].bm_response((tx,ty,tz),pol='x')[0]
#_bmy = aa[0].bm_response((tx,ty,tz),pol='y')[0]
bmxx = n.where(tz > 0, _bmx**2, 0)
#bmyy = n.where(tz > 0, _bmy**2, 0)
#bm_I = 0.5 * (bmxx + bmyy)
#bm_Q = 0.5 * (bmxx - bmyy)

#bl = n.array([100, 0, 0])
bl = aa.get_baseline(0,16,'r') * FQ

fig = p.figure(figsize=(5,5))
fig.subplots_adjust(left=.05, top=.95, bottom=.05, right=0.95)
m = Basemap(projection='ortho', lat_0=lat, lon_0=180, rsphere=1.)

def plot_hmap(dat, mode='log', mx=0, drng=3):
    d = proj_hmap(dat,'e')
    if mode == 'log': d = n.log10(n.abs(d).clip(10**(mx-drng),n.Inf))
    elif mode == 'real': d = d.real
    m.imshow(d, vmax=mx, vmin=mx-drng, origin='lower', interpolation='bicubic', alpha=.75, cmap='coolwarm')

fng = n.exp(2j*n.pi*(bl[0]*tx + bl[1]*ty + bl[2]*tz))
plot_hmap(fng, mode='real', mx=1, drng=2)
m.drawmapboundary(linewidth=2)
m.drawmeridians(n.arange(0, 360, 30), linewidth=2)
m.drawparallels(n.arange(-90,90,30), linewidth=2)

span = 30.
lons = n.arange(180-span/2-.005,180+span/2, .01)
for lt in n.arange(-90, 90, 10):
    lats = n.ones_like(lons) * lt
    x,y = m(lons,lats)
    m.plot(x,y, 'k', linewidth=1)
    if x[0]-x[1] == 0 and y[0]-y[1] == 0: continue
    print x[0]-x[1], y[0]-y[1]
    p.arrow(x[1],y[1],x[0]-x[1],y[0]-y[1], edgecolor='k', facecolor='k', overhang=.5, head_width=.04, length_includes_head=True)
p.show()
