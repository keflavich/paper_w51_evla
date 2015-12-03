import pyregion
import numpy as np
import pylab as pl
import paths
from paths import dpath,fpath,rpath
from spectral_cube import SpectralCube
import radio_beam
from astropy import coordinates
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy import log
import aplpy
from astropy.visualization import SqrtStretch,AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import common_constants
from common_constants import distance
pl.matplotlib.rc_file('pubfiguresrc')

files = {'natural': 'W51Ku_BD_h2co_v30to90_natural_contsub.image.fits',
         'briggs': 'W51Ku_BD_h2co_v30to90_briggs0_contsub.image.fits',
        }
levels = {'natural': [-0.4,-0.3,-0.2,-0.1, 0.010, 0.020, 0.030, 0.040],
          'briggs': [-0.28,-0.21,-0.14,-0.07, 0.007, 0.0105, 0.014],
         }

region = pyregion.open(rpath('w51e2zoom.reg'))

def set_tight_ticks(F):
    F.tick_labels.set_yformat('dd:mm:ss.ss')
    F.tick_labels.set_xformat('hh:mm:ss.ss')
    F.tick_labels.set_x_full_label_side('left')

for name, fn in files.iteritems():
    cube = SpectralCube.read(paths.dpath(fn))
    scube = (cube
             .subcube_from_ds9region(region)
             .with_spectral_unit(u.km/u.s, velocity_convention='radio')
             .spectral_slab(50*u.km/u.s, 63*u.km/u.s))
    beam = radio_beam.Beam.from_fits_header(cube.header)

    sigma = np.nanpercentile(scube.std(axis=0), 25)
    min = scube.min(axis=0)
    max = scube.max(axis=0)
    ok = (min.value < -5*sigma) | (max.value > 3*sigma)
    mask = (scube > sigma*scube.unit) | (scube < -5*sigma*scube.unit)

    m1 = scube.with_mask(mask).moment1(axis=0)
    hdu = m1.hdu
    hdu.data[~ok] = np.nan
    hdu.header.update(beam.to_header_keywords())

    m0 = scube.with_mask(mask).moment0(axis=0)

    pl.figure(1).clf()

    F = aplpy.FITSFigure(hdu)
    F.set_auto_refresh(False)
    set_tight_ticks(F)
    F.show_colorscale(cmap=pl.cm.RdYlBu_r, vmin=54, vmax=59)
    F.show_colorbar()
    F.add_beam()
    F.beam.set_facecolor('none')
    F.beam.set_hatch('//')
    F.refresh()
    F.save(fpath('velocity/w51e2zoom_{0}.png'.format(name)))

    F.show_regions(rpath('shi2010.reg'), layer='shi2010')
    for artist in F.get_layer('shi2010_txt').artistlist:
        artist.set_color('black')
    F.refresh()
    F.save(fpath('velocity/w51e2zoom_{0}_labeled.png'.format(name)))
    F.hide_layer('shi2010')
    F.hide_layer('shi2010_txt')
    F.show_regions(rpath('shi2010_notext.reg'), layer='shi2010_notext')
    for artist in F.get_layer('shi2010_notext_txt').artistlist:
        artist.set_color('black')
    F.refresh()
    F.save(fpath('velocity/w51e2zoom_{0}_marked.png'.format(name)))

    #F.show_contour(m0.hdu, levels=[-0.4,-0.3,-0.2,-0.1, 0.020, 0.040], colors=['k']*8)
    F.show_contour(m0.hdu, levels=[ L for L in levels[name] if L < 0], colors=['k']*8,
                   linewidths=[0.5]*8, alpha=0.5, linestyle='--')
    F.show_contour(m0.hdu, levels=[ L for L in levels[name] if L > 0], colors=['k']*8,
                   linewidths=[0.5]*8, alpha=0.5)
    F.refresh()
    F.save(fpath('velocity/w51e2zoom_{0}_contoured_marked.png'.format(name)))

    F.hide_layer('shi2010_notext')
    F.hide_layer('shi2010_notext_txt')
    F.refresh()
    F.save(fpath('velocity/w51e2zoom_{0}_contoured.png'.format(name)))
