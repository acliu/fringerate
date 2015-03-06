#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, capo as C, scipy
from mpl_toolkits.basemap import Basemap
from matplotlib import gridspec

#p.ion()
FQ = .150

cmap = p.get_cmap('jet')

aa = a.cal.get_aa('psa32_maxred', .1, FQ, 1)

print 'Plotting sky at %f GHz' % (FQ)
args = ['data/gsm/dOC-150MHz.fits']
print 'Reading %s' % args[0]
gsm = a.map.Map(fromfits=args[0])
print 'SCHEME:', gsm.scheme()
print 'NSIDE:', gsm.nside()
gsm.set_interpol(False)

DIM,RES = 400, .5
SIZE = DIM / RES
im = a.img.Img(DIM, res=RES)
h = a.healpix.HealpixMap(nside=gsm.nside())

def proj_hmap(dat, lst=0, dec=aa.lat, mode='e'):
    h.map = dat
    if mode == 'e' or mode == 'g': x,y,z = im.get_eq(lst, dec=dec, center=(SIZE/2,SIZE/2))
    elif mode == 't': x,y,z = im.get_top(center=(SIZE/2,SIZE/2))
    else: raise ValueError(mode)
    v = n.logical_not(x.mask)
    x,y,z = x.flatten(), y.flatten(), z.flatten()
    if mode == 'g':
        m = a.coord.convert_m('eq', 'ga', iepoch=a.ephem.J2000, oepoch=a.ephem.J2000)
        x,y,z = n.dot(m, n.array([x,y,z]))
    d = h[x,y,z]
    d.shape = (SIZE,SIZE)
    d = n.where(v, d, 0)
    return d

m = Basemap(projection='ortho', lat_0=90, lon_0=180, rsphere=1.)

def plot_hmap(dat, isys='e', lst=0., dec=aa.lat, mode='log', mx=None, drng=None):
    d = proj_hmap(dat, mode=isys, lst=lst, dec=dec)
    if mode == 'log': d = n.log10(n.abs(d))
    elif mode == 'dB':
        d = 10*n.log10(n.abs(d))
        mx,drng = 10*mx, 10*drng
    elif mode == 'real': d = d.real
    if mx is None: mx = d.max()
    if drng is None: drng = mx - d.min()
    print mx, drng
    d = d.clip(mx-drng,mx)
    im = m.imshow(d, vmax=mx, vmin=(mx-drng), origin='lower', interpolation='bicubic', cmap='jet')
    #m.drawmapboundary()
    #m.drawmeridians(n.arange(0,360,30))
    #m.drawparallels(n.arange(0,90,10))
    return im

def contour_hmap(dat, levels, isys='e', lst=0., dec=aa.lat, mode='log', colors='m', linestyles='-'):
    d = proj_hmap(dat, mode=isys, lst=lst, dec=dec)
    if mode == 'log': d = n.log10(n.abs(d))
    elif mode == 'dB':
        d = 10*n.log10(n.abs(d))
        mx,drng = 10*mx, 10*drng
    elif mode == 'real': d = d.real
    lons,lats,x,y = m.makegrid(SIZE,SIZE, returnxy=True)
    cs = m.contour(x,y, d, levels, linewidth=8, linestyles=linestyles, colors=colors)
    p.clabel(cs, inline=1, fontsize=10, fmt='%4.2f')#, manual=manual)
    #m.drawmapboundary()
    #m.drawmeridians(n.arange(0,360,30))
    #m.drawparallels(n.arange(0,90,10))
    return cs

xyz = h.px2crd(n.arange(h.npix()), ncrd=3)
top = n.dot(aa._eq2zen, xyz)
bl = aa.get_baseline(0,16,'r') * FQ
print 'Baseline:', bl
fng = C.frf_conv.mk_fng(bl, *xyz)
#plot_hmap(fng, mode='real'); p.show()

frf,bins,wgt,(cen,wid) = C.frf_conv.get_optimal_kernel_at_ref(aa, 0, (0,16))
bwfrs = C.frf_conv.get_beam_w_fr(aa, (0,16), ref_chan=0)
tbins,firs,frbins,frfs = C.frf_conv.get_fringe_rate_kernels(bwfrs,42.9,403)
wgt = scipy.interpolate.interp1d(frbins, frfs[0], kind='linear', bounds_error=False, fill_value=0)
fng_wgt = wgt(fng)

_bmx = aa[0].bm_response(top,pol='x')[0]
_bmy = aa[0].bm_response(top,pol='y')[0]
bm = 0.5 * (_bmx**2 + _bmy**2)
bm = n.where(top[-1] > 0, bm, 0)
fng_bm = fng_wgt * bm


fig = p.figure(figsize=(5,6.5))
gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
fig.subplots_adjust(left=.16, top=.95, bottom=.10, wspace=.05, hspace=.05, right=0.90)

p.subplot(gs[0])
plot_hmap(gsm.map.map, drng=2, isys='g', lst=1.5)
levels = [.5,.1,.01]
contour_hmap(bm, levels, mode='real', isys='e', linestyles='--', colors='k')#'#ff7f0e')
contour_hmap(fng_bm, levels, mode='real', isys='e', colors='r')

p.subplot(gs[1])
uv = a.miriad.UV('data/gsm/vis_simulation_1002.uv')
fqs = a.cal.get_freqs(uv['sdf'],uv['sfreq'],uv['nchan'])
uvL = a.miriad.UV('data/gsm/vis_simulation_1002.uvL')
a.scripting.uv_selector(uv, '0_16', 'xx')
a.scripting.uv_selector(uvL, '0_16', 'xx')
ds, dLs = [],[]
for pre,d,f in uv.all(raw=True):
    if n.abs(uv['lst'] - 1.5) < .01: ds.append(d)
for pre,d,f in uvL.all(raw=True):
    if n.abs(uvL['lst'] - 1.5) < .01: dLs.append(d)
ds,dLs = n.array(ds), n.array(dLs)
print ds.shape, fqs.shape
ds = n.sqrt(n.average(n.abs(ds)**2, axis=0))
dLs = n.sqrt(n.average(n.abs(dLs)**2, axis=0))
print ds.shape
p.semilogy(1e3*fqs, ds, 'k')
p.semilogy(1e3*fqs, dLs, 'r')
p.xlim(1e3*fqs[0], 1e3*fqs[-1])
p.xlabel('Frequency [MHz]')
p.ylabel('RMS Visiblity [K]')
p.grid()

p.show()

