import os
import numpy as np
import pylab as pl
from paths import dpath,fpath,rpath
from spectral_cube import SpectralCube
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
import aplpy
import paths
from common_constants import distance
import FITS_tools

pl.close(1)
figure = pl.figure(1)
figure.clf()

hdu = fits.open(paths.dpath('naco_Kband_W51.fits'))

F = aplpy.FITSFigure(hdu,convention='calabretta',figure=figure)
F.set_auto_refresh(False)
#F = aplpy.FITSFigure(dpath+'W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.rot45.fits',convention='calabretta',figure=figure)
F.tick_labels.set_xformat('dd.ddd')
F.tick_labels.set_yformat('dd.ddd')
F.tick_labels.set_font(size=20)
F.axis_labels.set_font(size=20)
F.show_grayscale(stretch='log',vmin=-1,vmid=-1.5,vmax=1e2)
#e1 = coordinates.ICRS(290.93263,14.50745,unit=('deg','deg'))
#F.recenter(e1.ra.value,e1.dec.value,width=1/60.,height=1/60.)
#F.recenter(290.92633,14.514769,radius=1.4/60.)

F.add_scalebar(length=((0.5 * u.pc)/distance*u.radian).to(u.degree).value)

F.recenter(290.91669,14.518151,radius=8./3600.)
F.recenter(290.91766,14.518604,width=20./3600.,height=10./3600.)
F.scalebar.set_length(((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)


h77a = SpectralCube.read(paths.dpath('W51north_H77_Outflow_cutout.fits'))
h77a_outflow = h77a.spectral_slab(-16*u.km/u.s, -60*u.km/u.s).sum(axis=0)
c = (0,0.1,0.9)
h77acolors = [c[:3] + (1,) for x in (0.1,0.2,0.3,0.4,0.5,0.6,0.7)]
h77alevels = np.arange(0.02,0.05,0.005)
F.show_contour(h77a_outflow.hdu, levels=h77alevels, colors=h77acolors,
               filled=False, layer='h77a')

cutout_coords = {'xlo':290.91974*u.deg, 'xhi':290.90926*u.deg, 'ylo':14.51492*u.deg, 'yhi':14.523815*u.deg}
cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio').subcube(**cutout_coords)

vr = [56,60]
core = cube.spectral_slab(vr[0]*u.km/u.s, vr[1]*u.km/u.s).sum(axis=0)
c = (0.9,0.1,0.0)
corecolors = [c[:3] + (1,) for x in (0.1,0.2,0.3,0.4,0.5,0.6,0.7)]
corelevels = np.arange(0.02,0.05,0.005)
F.show_contour(core.hdu, levels=corelevels, colors=corecolors,
               filled=False, layer='core')

F.show_regions(paths.rpath('w51_sinfoni_pointings.reg'))
F.save(os.path.expanduser('~/proposals/vlt/p96/w51/SINFONI_pointings_on_NACO.png'))

pl.draw(); pl.show()
