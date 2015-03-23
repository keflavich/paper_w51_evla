import numpy as np
import pylab as pl
from paths import dpath,fpath,rpath
from spectral_cube import SpectralCube
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
import aplpy
import pylab

import pyregion
import gaussfitter

cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')

r = pyregion.open(rpath('W51_22_emission.reg'))
scube = cube.subcube_from_ds9region(pyregion.ShapeList([x for x in r if 'text' in x.attr[1] and x.attr[1]['text']=='W51NorthCore']))


vr = [56,61]
gaussfits = [gaussfitter.gaussfit(
    np.nan_to_num((scube[scube.closest_spectral_channel(v*u.km/u.s),:,:]/scube.max()).value),
    return_error=True )  for v in np.arange(vr[0],vr[1]+0.5,0.5)]

centroids = np.array([gf[0][2:4] for gf in gaussfits])
ecentroids = np.array([gf[1][2:4] for gf in gaussfits])

fig = pl.figure(1)
fig.clf()
ax = fig.gca()
sc = ax.scatter(centroids[:,0], centroids[:,1], marker='o', c=np.arange(vr[0],vr[1]+0.5,0.5), s=100)

# sc doesn't have facecolors until it is drawn
pl.draw()
itr = zip(centroids[:,0], centroids[:,1], ecentroids[:,0], ecentroids[:,1], sc.get_facecolors())

print itr
print

for x,y,ex,ey,c in itr:
    print x,y,ex,ey,c
    ax.errorbar(x,y,xerr=ex,yerr=ey, linestyle='none', zorder=-1, color=c)

fig.colorbar(sc)
fig.savefig(fpath('w51north_core_gaussfits.png'))
pl.draw()
pl.show()

