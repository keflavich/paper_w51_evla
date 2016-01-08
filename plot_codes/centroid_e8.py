"""
results:
natural
Integrated gaussfit:  [  1.36949123e-02   1.93799694e-02   6.11740274e+00   5.22358496e+00
   2.57607937e+00   3.27451343e+00   1.55816943e+02]
Center position: 19h23m43.9059s +14d30m28.241s
cube max: 0.007375472225248814 Jy
beam: 2.624705823202675e-11 sr
nbeams: 1.8985191167945263
Integrated intensity: 18.396621239362144 km mJy / s
Integrated brightness: 108.68478434233081 K km / s
Integrated average brightness: 57.24713719282164 K km / s
radius: 0.5808754255642093 arcsec = 0.015235444715361425 pc

briggs0
Integrated gaussfit:  [  4.59430179e-03   1.05029892e-02   1.40664908e+01   1.02741558e+01
   2.02167759e+00   3.74588596e+00   1.67305242e+02]
Center position: 19h23m43.8933s +14d30m28.2237s
cube max: 0.0036180811002850533 Jy
beam: 4.7859023011910344e-12 sr
nbeams: 2.3368574507226474
Integrated intensity: 12.271994330020272 km mJy / s
Integrated brightness: 397.61488396662423 K km / s
Integrated average brightness: 170.14939608047817 K km / s
radius: 0.275190364769809 arcsec = 0.007217808507871056 pc
"""

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

for name, filename in (('natural','W51Ku_BD_h2co_v30to90_natural_contsub.image.fits'),
                       ('briggs0','W51Ku_BD_h2co_v30to90_briggs0_contsub.image.fits')):

    cube_ = SpectralCube.read(dpath(filename))
    cube = cube_.with_spectral_unit(u.km/u.s, velocity_convention='radio')
    beam = radio_beam.Beam.from_fits_header(cube.header)

    r = pyregion.open(rpath('W51_22_emission.reg'))
    r_subset = pyregion.ShapeList([x for x in r if 'text' in x.attr[1] and
                                   x.attr[1]['text']=='e8mol_ext'])
    scube = cube.subcube_from_ds9region(r_subset)
    pixscale = scube.wcs.pixel_scale_matrix[1,1]*u.deg
    pixscale_as = pixscale.to(u.arcsec)

    cubemax = scube.max() # normalize everything
    noise = scube.spectral_slab(30*u.km/u.s, 50*u.km/u.s).std(axis=0)
    noise[np.isnan(noise.value)] = 1e10*noise.unit
    noise_norm = (noise / cubemax).value

    dv = cube.spectral_axis.diff().mean()
    vr = [57,64]
    gaussfits = [gaussfitter.gaussfit(
        np.nan_to_num((scube[scube.closest_spectral_channel(v*u.km/u.s),:,:]/cubemax).value),
        err=noise_norm,
        return_error=True)
        for v in ProgressBar(np.arange(vr[0],vr[1]+dv.value,dv.value))]

    integrated_image = (scube.spectral_slab(vr[0]*u.km/u.s,
                                            vr[1]*u.km/u.s).moment0(axis=0))
    gaussfit_total,fitimage = gaussfitter.gaussfit(np.nan_to_num(integrated_image.value),
                                                   err=noise.value,
                                                   returnfitimage=True
                                         )
    # sanity check: make sure the fit is OK (it is)
    pl.figure(2).clf()
    pl.imshow(integrated_image.value, cmap=pl.cm.bone_r)
    cb = pl.colorbar()
    pl.contour(fitimage, cmap=pl.cm.spectral)
    pl.savefig(fpath('sanitycheck_w51e8_core_gaussfits_{0}.png'.format(name)))

    pl.figure(3).clf()
    pl.imshow((u.Quantity(integrated_image)/(u.km/u.s))
              .to(u.K, u.brightness_temperature(beam, 14.488*u.GHz)).value,
              cmap=pl.cm.bone_r)
    cb = pl.colorbar()
    pl.contour(fitimage, cmap=pl.cm.spectral)
    pl.savefig(fpath('sanitycheck_K_w51e8_core_gaussfits_{0}.png'.format(name)))

    print(name)
    print("Integrated gaussfit: ", gaussfit_total)
    centerx, centery = gaussfit_total[2:4]
    celwcs = scube.wcs.sub([wcs.WCSSUB_CELESTIAL])
    cx_w, cy_w = celwcs.wcs_pix2world([[centerx, centery]], 0)[0]
    center_coord = coordinates.SkyCoord(cx_w, cy_w, unit=('deg','deg'), frame='fk5')
    print("Center position: {0}".format(center_coord.to_string('hmsdms')))
    integ_intens = (1.*u.Jy/beam * gaussfit_total[1] *
                    (gaussfit_total[4]*pixscale_as *
                     gaussfit_total[5]*pixscale_as * 2 * np.pi))
    nbeams = ((gaussfit_total[4]*pixscale_as * gaussfit_total[5]*pixscale_as * 2 * np.pi) / beam).decompose()
    print("cube max: {0}".format(cubemax))
    print("beam: {0}".format(beam))
    print("nbeams: {0}".format(nbeams))
    print("Integrated intensity: {0}".format(integ_intens.to(u.mJy)*dv))
    kkms = integ_intens.to(u.K, u.brightness_temperature(beam, 14.488*u.GHz))*dv
    print("Integrated brightness: {0}".format(kkms))
    mean_kkms = integ_intens.to(u.K, u.brightness_temperature(beam,
                                                              14.488*u.GHz))*dv/nbeams
    print("Integrated average brightness: {0}".format(mean_kkms))
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
    fig.savefig(fpath('w51e8_core_gaussfits_{0}.png'.format(name)))
    pl.draw()
    pl.show()

    print()
