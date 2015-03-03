#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, re
import glob, os
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_bm_hmap.py [options] *.hmap')
opts,args = o.parse_args(sys.argv[1:])

srcfiles = sys.argv[1:]

def beam_area(dat):
    return n.sum(dat) * 4 * n.pi / dat.size
def beamsq_area(dat):
    return n.sum(dat**2) * 4 * n.pi / dat.size

FQ = .151
DIM = 400
RES = .5
SIZE = DIM / RES

aa = a.cal.get_aa('psa898_v003', n.array([FQ]))

h = a.healpix.HealpixMap(nside=64)
im = a.img.Img(DIM, res=RES)

def proj_hmap(dat, mode='t'):
    h.map = dat
    if mode == 'e': x,y,z = im.get_eq(dec=aa.lat, center=(SIZE/2,SIZE/2))
    elif mode == 't': x,y,z = im.get_top(center=(SIZE/2,SIZE/2))
    else: raise ValueError(mode)
    v = n.logical_not(x.mask)
    x,y,z = x.flatten(), y.flatten(), z.flatten()
    #h.set_interpol(True)
    d = h[x,y,z]
    d.shape = (SIZE,SIZE)
    d = n.where(v, d, 0)
    return d

m = Basemap(projection='ortho', lat_0=90, lon_0=180, rsphere=1.)

def plot_hmap(dat, mode='log', mx=0, drng=3):
    d = proj_hmap(dat,'t')
    d /= d.max() # control for any flux scale change in fringe rate filter
    if mode == 'log': d = n.log10(n.abs(d).clip(10**(mx-drng),n.Inf))
    elif mode == 'real': d = d.real
    #m.imshow(d, vmax=mx, vmin=mx-drng, origin='lower', interpolation='bicubic', cmap='jet')
    m.imshow(10*d, vmax=10*mx, vmin=10*(mx-drng), origin='lower', interpolation='bicubic', cmap='jet')
    m.drawmapboundary()
    m.drawmeridians(n.arange(0,360,30))
    m.drawparallels(n.arange(0,90,10))

gps = {
    1:glob.glob('data/d*/*_sim.npz'),
    2:glob.glob('data/flat_beam/d*/*_fng.npz'),
    3:glob.glob('data/d*/*_fng.npz'),
}

p.figure(figsize=(10,3))
for gp in gps:
    hmaps = [a.map.Map(nside=nside) for nside in [64,32,16,8]]

    for filename in gps[gp]:
        print '    Reading', filename
        npz = n.load(filename)
        fq = npz['freq']
        print '    Freq:', fq
        sdat = n.abs(npz['spec'])
        valid = n.where(sdat == 0, 0, 1)
        sdat = sdat.compress(valid)
        sx = npz['x'].compress(valid)
        sy = npz['y'].compress(valid)
        sz = npz['z'].compress(valid)
        #p.plot(sz,sdat); p.show()
        swgt = n.ones_like(sdat)
        for hmap in hmaps: hmap.add((sx,sy,sz), swgt, sdat)

    for hmap in hmaps: hmap.reset_wgt()
    d,w = n.zeros(hmaps[0].npix()), n.zeros(hmaps[0].npix())
    mx,my,mz = hmaps[0].px2crd(n.arange(hmaps[0].npix()))
    for hmap in hmaps:
        #hmap.set_interpol(True)
        wgt = hmap.nside()
        w += hmap.wgt[mx,my,mz] * wgt
        d += hmap.map[mx,my,mz] * wgt
    data  = n.where(w > 0, d/w, 0)
    hmap = hmaps[0]
    hmap.map.map,hmap.wgt.map = d,w
    p.subplot(1,3,gp)
    plot_hmap(hmap, mx=0, drng=2)
    p.subplots_adjust(left=.05,right=.87,top=.9,bottom=.1, wspace=.05)
fig = p.gcf()
cbar_ax1 = fig.add_axes([.9,.15,.02,.7])
p.colorbar(ticks=[0,-3,-6,-10,-13,-16,-20], cax=cbar_ax1)#, shrink=.5)
p.xlabel('dB', fontsize=14)

#Omega_P, Omega_PP = beam_area(data), beamsq_area(data)
#print 'Omega_P:', Omega_P
#print 'Omega_PP:', Omega_PP
#print 'Omega_eff:', Omega_P**2/Omega_PP
#bl_ew_len = 100. # ns
#fng_mx = (fq * bl_ew_len * 2*n.pi / a.const.sidereal_day)
#fng_mn = -0.5 * fng_mx
#fng_bins = n.linspace(fng_mn, fng_mx, 1000)
#fng_wgts = fng_wgt(fng_bins, fq, bl_ew_len)
#t_eff = 1. / n.average(fng_wgts**2)
#signal_attenuation = Omega_PP / 0.323
#print 't_fng/t_sky:', t_eff
#print 'signal attenuation:', signal_attenuation
#print 'Delta^2_N,fng / Delta^2_N,sky single mode:', 1. / signal_attenuation / t_eff
#print 'Delta^2_N,fng / Delta^2_N,sky avg modes:', 1. / signal_attenuation / n.sqrt(t_eff)


p.show()
