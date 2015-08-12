import numpy as np
import pylab as pl
from paths import dpath,fpath,rpath
from spectral_cube import SpectralCube
from astropy.utils.console import ProgressBar
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
from astropy import wcs
import aplpy
import pylab

import pyregion
import gaussfitter

cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')

r = pyregion.open(rpath('W51_22_emission.reg'))
scube = cube.subcube_from_ds9region(pyregion.ShapeList([x for x in r if 'text' in x.attr[1] and x.attr[1]['text']=='NorthCore']))
pixscale = scube.wcs.pixel_scale_matrix[1,1]*u.deg
pixscale_as = pixscale.to(u.arcsec)

cubemax = scube.max() # normalize everything
noise = (scube.spectral_slab(30*u.km/u.s, 50*u.km/u.s).std(axis=0) / cubemax).value
noise[np.isnan(noise)] = 1e10

vr = [56,63.5]
gaussfits = [gaussfitter.gaussfit(
    np.nan_to_num((scube[scube.closest_spectral_channel(v*u.km/u.s),:,:]/cubemax).value),
    err=noise,
    return_error=True )  for v in ProgressBar(np.arange(vr[0],vr[1]+0.5,0.5))]

gaussfit_total = gaussfitter.gaussfit(np.nan_to_num((scube.spectral_slab(vr[0]*u.km/u.s, vr[1]*u.km/u.s).mean(axis=0)/cubemax).value))
centerx, centery = gaussfit_total[2:4]
center_coord = coordinates.SkyCoord(*scube.wcs.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([[centerx, centery]], 0)[0], unit=('deg','deg'), frame='fk5')
print("Center position: {0}".format(center_coord.to_string('hmsdms')))

centroids = np.array([gf[0][2:4] for gf in gaussfits])
ecentroids = np.array([gf[1][2:4] for gf in gaussfits])

fig = pl.figure(1)
fig.clf()
ax = fig.gca()
sc = ax.scatter((centroids[:,0]-centerx)*pixscale_as.value,
                (centroids[:,1]-centery)*pixscale_as.value, marker='o',
                c=np.arange(vr[0],vr[1]+0.5,0.5), s=100, edgecolor='none',
                cmap=pl.cm.spectral)

# sc doesn't have facecolors until it is drawn
pl.draw()
itr = zip(centroids[:,0], centroids[:,1], ecentroids[:,0], ecentroids[:,1], sc.get_facecolors())

print itr
print

for x,y,ex,ey,c in itr:
    print x,y,ex,ey,c
    ax.errorbar((x-centerx)*pixscale_as.value,
                (y-centery)*pixscale_as.value,
                xerr=ex*pixscale_as.value,
                yerr=ey*pixscale_as.value, linestyle='none', zorder=-1, color='k')
pl.xlabel("Offset from centroid (arcsec)")
pl.ylabel("Offset from centroid (arcsec)")

cb = fig.colorbar(sc)
cb.set_label("$V_{LSR}$")
fig.savefig(fpath('w51north_core_gaussfits.png'))
pl.draw()
pl.show()

