#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C, scipy
from mpl_toolkits.basemap import Basemap

FQ = .151
DIM = 400
RES = .5
SIZE = DIM / RES

aa = a.cal.get_aa('psa898_v003', n.array([FQ]))

h = a.healpix.HealpixMap(nside=64)
im = a.img.Img(DIM, res=RES)

def proj_hmap(dat, mode='e'):
    h.map = dat
    if mode == 'e': x,y,z = im.get_eq(dec=aa.lat, center=(SIZE/2,SIZE/2))
    elif mode == 't': x,y,z = im.get_top(center=(SIZE/2,SIZE/2))
    else: raise ValueError(mode)
    v = n.logical_not(x.mask)
    x,y,z = x.flatten(), y.flatten(), z.flatten()
    d = h[x,y,z]
    d.shape = (SIZE,SIZE)
    d = n.where(v, d, 0)
    return d

#m = Basemap(projection='ortho', lat_0=aa.lat, lon_0=180, rsphere=1.)
m = Basemap(projection='ortho', lat_0=90, lon_0=180, rsphere=1.)

def plot_hmap(dat, mode='log', mx=0, drng=3):
    d = proj_hmap(dat,'e')
    #d /= d.max()
    if mode == 'log': d = n.log10(n.abs(d).clip(10**(mx-drng),n.Inf))
    elif mode == 'dB':
        d = 10*n.log10(n.abs(d).clip(10**(mx-drng),n.Inf))
        mx,drng = 10*mx, 10*drng
    elif mode == 'real': d = d.real
    im = m.imshow(d, vmax=mx, vmin=(mx-drng), origin='lower', interpolation='bicubic', cmap='jet')
    m.drawmapboundary()
    m.drawmeridians(n.arange(0,360,30))
    m.drawparallels(n.arange(0,90,10))
    return im

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
#p.plot(bins, hXX); p.plot(bins, hYY); p.show()
#p.plot(bins, n.where(hXX>0,hYY/hXX,0)); p.show()
#bm_pol = n.sqrt(polmatch(fng))
bm_pol = polmatch(fng)
bmXXm = bmXX 
bmYYm = bmYY * bm_pol

hpol, bin_edges = n.histogram(fng, bins=bin_edges, weights=(bmXXm-bmYYm)**2)
htot, bin_edges = n.histogram(fng, bins=bin_edges, weights=(bmXXm+bmYYm)**2)
bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
#match = n.where(hpol > 0, n.sqrt(htot/hpol), 0); match /= match.max() # XXX sqrt here or no?
match = n.where(hpol > 0, htot/hpol, 0); match /= match.max() # XXX sqrt here or no?
wgt_pol = scipy.interpolate.interp1d(bins, match, kind='linear', bounds_error=False, fill_value=0)
tbins,firs,frbins,frfs= C.frf_conv.get_fringe_rate_kernels(n.array([wgt_pol]),42.9,403)
frfs[0] /= frfs[0].max()
wgt_pol = scipy.interpolate.interp1d(frbins, frfs[0], kind='linear')


frf,bins,wgt,(cen,wid) = C.frf_conv.get_optimal_kernel_at_ref(aa, 0, (0,16))

bwfrs = C.frf_conv.get_beam_w_fr(aa, (0,16), ref_chan=0) # XXX check our ref_chan
print bwfrs
tbins,firs,frbins,frfs = C.frf_conv.get_fringe_rate_kernels(bwfrs,42.9,403)
wgt = scipy.interpolate.interp1d(frbins, frfs[0], kind='linear')
#p.plot(bins,wgt(bins))
#p.plot(bins,wgt_pol(bins))
#p.plot(bins,wgt(bins)*wgt_pol(bins)); p.show()
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

frf_fngpol = wgt_pol(frbins) * wgt(frbins)
skypass = n.where(n.abs(frbins) < .0012, 1., 0)
crit_noise_lev = n.sum(skypass**2) / frbins.size
fng_noise_lev = n.sum(n.abs(frf_fngpol)**2)/ frbins.size
print 'Eff Int Time:', 42.9/crit_noise_lev, 42.9/fng_noise_lev
print 'Noise Redux:', fng_noise_lev / crit_noise_lev
print 'Sensitivity:', n.sqrt(crit_noise_lev / fng_noise_lev) * n.sum(bmIm_fngpol**2)/n.sum(bmI**2)

print beam_area(bmI), beam_area(bmIm), beam_area(bmIm_fng), beam_area(bmIm_fngpol)
print beam_area(bmI**2), beam_area(bmIm**2), beam_area(bmIm_fng**2), beam_area(bmIm_fngpol**2)

print 'XXX', beam_area(bmI**2), beam_area((bmI*fng_wgt)**2)

print 'Unweighted:', n.sqrt(n.sum(bmQ**2)/n.sum(bmI**2))
print 'Matched:   ', n.sqrt(n.sum(bmQm**2)/n.sum(bmIm**2))
print 'Optimal SN:', n.sqrt(n.sum(bmQm_fng**2)/n.sum(bmIm_fng**2))
print 'Opt pol+SN:', n.sqrt(n.sum(bmQm_fngpol**2)/n.sum(bmIm_fngpol**2))

print 'Unweighted:', (n.sum(bmQ**2)/n.sum(bmI**2))
print 'Matched:   ', (n.sum(bmQm**2)/n.sum(bmIm**2))
print 'Optimal SN:', (n.sum(bmQm_fng**2)/n.sum(bmIm_fng**2))
print 'Opt pol+SN:', (n.sum(bmQm_fngpol**2)/n.sum(bmIm_fngpol**2))

print 'Whole P beam area:', beam_area(bmI)
print 'Matched P beam area:', beam_area(bmIm)
print 'Optimal P beam area:', beam_area(bmIm_fng)
print 'Opt pol P beam area:', beam_area(bmIm_fngpol)

print 'Whole PP beam area:', beam_area(bmI**2)
print 'Matched PP beam area:', beam_area(bmIm**2)
print 'Optimal PP beam area:', beam_area(bmIm_fng**2)
print 'Opt pol PP beam area:', beam_area(bmIm_fngpol**2)

p.figure(figsize=(12,6))
kw = {'mode':'dB','mx':0,'drng':2}
p.subplot(241); plot_hmap(bmI, **kw)
p.subplot(242); plot_hmap(bmIm, **kw)
p.subplot(243); plot_hmap(bmIm_fng, **kw)
p.subplot(244); im1 = plot_hmap(bmIm_fngpol, **kw)

kw = {'mode':'real','mx':.05,'drng':.1}
p.subplot(245); plot_hmap(bmQ, **kw)
p.subplot(246); plot_hmap(bmQm, **kw)
p.subplot(247); plot_hmap(bmQm_fng, **kw)
p.subplot(248); im2 = plot_hmap(bmQm_fngpol, **kw)

p.subplots_adjust(left=.05,right=.87,top=.9,bottom=.1, wspace=.05, hspace=.05)
fig = p.gcf()
cbar_ax1 = fig.add_axes([.9,.56,.02,.3])
p.colorbar(im1, ticks=[0,-3,-6,-10,-13,-16,-20], cax=cbar_ax1)#, shrink=.5)
p.xlabel(r'$\rm dB$', fontsize=14)
cbar_ax2 = fig.add_axes([.9,.13,.02,.3])
p.colorbar(im2, cax=cbar_ax2)#, shrink=.5)
p.xlabel(r'$({\rm XX}-{\rm YY})/2$', fontsize=14)
p.show()

#def nos(dat, t):
#    pwr_bm = beam_area(bm)
#    pwr2_bm = beam_area(dat**2)
#    return pwr_bm**2 / pwr2_bm / t

#print beam_area(bm), beam_area(bm**2), nos(bm, 1.), nos(bm,1.) / n.sqrt(NBINS)

#xyz = (xyz[1],xyz[0],xyz[2])
#bl = n.array([100, 0, 0])
#fng = mk_fng(bl, *xyz)
print fng.max()


skypass = n.where(n.abs(bins)<0.5/42.9, 1, 0)
crosstalk = n.ones_like(frbins)
#crosstalk[14:] = 0 # 600s / 42.9s = 14
#crosstalk[3*14:] = 0 # 600s / 42.9s = 14
#crosstalk[:3*14] *= a.dsp.gen_window(3*14, 'blackman-harris')
#crosstalk[3*84:] = 0 # 3600s / 42.9s = 84
#crosstalk[:3*84] *= a.dsp.gen_window(3*84, 'blackman-harris')
crosstalk = n.fft.fftshift(n.abs(n.fft.ifft(crosstalk)))
crosstalk /= crosstalk.max()
print frfs[0].max()
frfs[0] /= frfs[0].max() # XXX shouldn't have to do this
#p.plot(frbins,1-crosstalk)
#p.plot(bins,frf)
#p.plot(frbins,frfs[0])
#frfs[0] *= (1-crosstalk)
#p.show()
#print fng.max(), fng.min()
noise_ratio = n.sum(n.abs(frfs[0])**2)/n.sum(n.abs(skypass)**2)
print noise_ratio, 42.9 / noise_ratio, n.sqrt(noise_ratio)


#hist, bin_edges = n.histogram(fng, bins=bin_edges, weights=bmI**2)
#hist, bin_edges = n.histogram(fng, range=(-.01,.01), bins=NBINS, weights=bm**2) # or bm**2?
#f = 0
#for cnt,b1 in enumerate(bin_edges[:-1]):
#    b2 = bin_edges[cnt+1]
#    if cnt % 2 == 0: f = n.where(n.logical_and(fng >= b1, fng < b2), 1., f)


#
#p.subplot(234)
##hist, bins = n.histogram(fng, range=(-.01,.01), bins=NBINS, weights=bm) # or bm**2?
#bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
#hist1 = hist/n.sum(hist)
#print 'Hist:', hist1.max()
##p.plot(bins, hist1/hist1.max())
#p.plot(bins, hist1/hist1.max())
#p.plot(bins, n.sqrt(hist1/hist1.max()))
#
#noise_sum, noise_wgt = 0, 0
#wgts = []
#f = 0
#for cnt,b1 in enumerate(bin_edges[:-1]):
#    b2 = bin_edges[cnt+1]
#    fng_wgt = n.where(n.logical_and(fng >= b1, fng < b2), 1., 0)
#    fng_amp = bm * fng_wgt
#    _n = nos(fng_amp,NBINS)
#    w = 1. / _n**2
#    wgts.append(w)
#    if w > 0:
#        #print cnt, b1, _n, w
#        noise_sum += (_n*w)**2
#        noise_wgt += w
#    f = n.where(n.logical_and(fng >= b1, fng < b2), w, f)
#wgts = n.array(wgts)
#
#p.plot(bins, wgts/wgts.max())
#
#fng_amp = bm * f
#print fng_amp.max()
#fng_amp /= fng_amp.max()
#print n.sqrt(noise_sum) / noise_wgt
#
#print 'Improvement:', (n.sqrt(noise_sum)/noise_wgt) / (nos(bm,1.) / n.sqrt(NBINS))
#    
#p.subplot(235)
#p.imshow(proj_hmap(f,'e'), origin='lower')
#
#p.subplot(236)
#p.imshow(proj_hmap(fng_amp,'e'), origin='lower')

p.show()
