import numpy as np
import pylab as pl
from paths import dpath,fpath,rpath
from spectral_cube import SpectralCube
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
import aplpy
from wcsaxes import WCS

distance = 5.41*u.kpc
pl.figure(1).clf()

uwish_hdu = fits.open(dpath('UKIDSS/113_215_30_460_1_hdu1.fits'))[1]

F = aplpy.FITSFigure(uwish_hdu, figure=pl.figure(1))
F.show_grayscale(vmax=4000, invert=True)
F.tick_labels.set_xformat('dd.ddddd')
F.tick_labels.set_yformat('dd.ddddd')
#F.ticks.set_yspacing(0.1)
#F.ticks.set_xspacing(0.1)
F.recenter(290.91557, 14.524966, radius=0.1/60.)

F.add_scalebar(length=((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('black')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)

hdu = fits.open(dpath('W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'))
for ii in [3,4]:
    for kk in ['CRVAL','CTYPE','CDELT','CRPIX','CUNIT','NAXIS']:
        k = kk+str(ii)
        if k in hdu[0].header:
            del hdu[0].header[k]
hdu[0].header['NAXIS'] = 2
hdu[0].data = hdu[0].data.squeeze()

#hdu = FITS_tools.hcongrid.hcongrid_hdu(hdu, fits.Header.fromtextfile('d4_cutout.hdr'))

F.show_contour(hdu, levels=[5e-5, 1e-4, 2e-4, 4e-4, 8e-4,],
               colors=['r']*5)

F.save(fpath('d4_h2_overlay.pdf'))
F.save(fpath('d4_h2_overlay.png'))


# still no ticks.
#mywcs = WCS(uwish_hdu.header)
#fig = pl.figure(2)
#fig.clf()
#ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=mywcs)
#ax.imshow(uwish_hdu.data, vmin=3285, vmax=4000, cmap=pl.cm.gray_r,
#          origin='lower')
#lon = ax.coords['ra']
#lon.set_major_formatter('hh:mm:ss.ssss')
#lat = ax.coords['dec']
#lat.set_major_formatter('dd:mm:ss.ssss')
#lon.set_ticks(number=3)
#lat.set_ticks(number=3)
#fig.savefig(fpath('d4_h2_overlay_wcsaxes.png'))
#
#ax.contour(hdu[0].data, levels=[5e-5, 1e-4, 2e-4, 4e-4, 8e-4,], colors=['r']*5,
#           transform=ax.get_transform(WCS(hdu[0].header)))
#fig.savefig(fpath('d4_h2_overlay_wcsaxes.png'))
