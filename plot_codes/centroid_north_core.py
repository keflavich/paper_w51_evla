import numpy as np
import pylab as pl
import common_constants
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
import radio_beam

cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
beam = radio_beam.Beam.from_fits_header(cube.header)

r = pyregion.open(rpath('W51_22_emission.reg'))
scube = cube.subcube_from_ds9region(pyregion.ShapeList([x for x in r if 'text' in x.attr[1] and x.attr[1]['text']=='NorthCore']))
pixscale = scube.wcs.pixel_scale_matrix[1,1]*u.deg
pixscale_as = pixscale.to(u.arcsec)

cubemax = scube.max() # normalize everything
noise = scube.spectral_slab(30*u.km/u.s, 50*u.km/u.s).std(axis=0)
noise[np.isnan(noise)] = 1e10*noise.unit
noise_norm = (noise / cubemax).value

dv = cube.spectral_axis.diff().mean()
vr = [56,63.5]
gaussfits = [gaussfitter.gaussfit(
    np.nan_to_num((scube[scube.closest_spectral_channel(v*u.km/u.s),:,:]/cubemax).value),
    err=noise_norm,
    return_error=True )  for v in ProgressBar(np.arange(vr[0],vr[1]+dv.value,dv.value))]

integrated_image = (scube.spectral_slab(vr[0]*u.km/u.s, vr[1]*u.km/u.s).sum(axis=0))
gaussfit_total,fitimage = gaussfitter.gaussfit(np.nan_to_num(integrated_image.value),
                                               err=noise.value,
                                               returnfitimage=True
                                     )
# sanity check: make sure the fit is OK (it is)
pl.figure(2).clf()
pl.imshow(integrated_image.value, cmap=pl.cm.bone_r)
pl.contour(fitimage, cmap=pl.cm.spectral)

print("Integrated gaussfit: ", gaussfit_total)
centerx, centery = gaussfit_total[2:4]
center_coord = coordinates.SkyCoord(*scube.wcs.sub([wcs.WCSSUB_CELESTIAL]).wcs_pix2world([[centerx, centery]], 0)[0], unit=('deg','deg'), frame='fk5')
print("Center position: {0}".format(center_coord.to_string('hmsdms')))
integ_intens = 1.*u.Jy/beam * gaussfit_total[1] * (gaussfit_total[4]*pixscale_as * gaussfit_total[5]*pixscale_as * 2 * np.pi)
nbeams = ((gaussfit_total[4]*pixscale_as * gaussfit_total[5]*pixscale_as * 2 * np.pi) / beam).decompose()
print("Integrated intensity: {0}".format(integ_intens.to(u.mJy)*dv))
print("Integrated brightness: {0}".format(integ_intens.to(u.K, u.brightness_temperature(beam, 14.488*u.GHz))*dv))
print("Integrated average brightness: {0}".format(integ_intens.to(u.K, u.brightness_temperature(beam, 14.488*u.GHz))*dv/nbeams))
radius = (gaussfit_total[4]*pixscale_as * gaussfit_total[5]*pixscale_as)**0.5
print("radius: {0} = {1}".format(radius,
                                 (radius*common_constants.distance).to(u.pc,
                                                                       u.dimensionless_angles())))

centroids = np.array([gf[0][2:4] for gf in gaussfits])
ecentroids = np.array([gf[1][2:4] for gf in gaussfits])

fig = pl.figure(1)
fig.clf()
ax = fig.gca()
sc = ax.scatter((centroids[:,0]-centerx)*pixscale_as.value,
                (centroids[:,1]-centery)*pixscale_as.value, marker='o',
                c=np.arange(vr[0],vr[1]+dv.value,dv.value), s=100, edgecolor='none',
                cmap=pl.cm.spectral)

# sc doesn't have facecolors until it is drawn
pl.draw()
itr = zip(centroids[:,0], centroids[:,1], ecentroids[:,0], ecentroids[:,1], sc.get_facecolors())

#print itr
#print

for x,y,ex,ey,c in itr:
    #print x,y,ex,ey,c
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

