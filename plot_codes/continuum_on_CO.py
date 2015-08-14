import paths
import pyregion
import numpy as np
import pylab as pl
from paths import dpath,fpath,rpath
from spectral_cube import SpectralCube
from astropy import coordinates
from astropy.utils.console import ProgressBar
from astropy import units as u
from astropy.io import fits
from astropy import wcs
from astropy import log
import aplpy
from astropy.visualization import SqrtStretch,AsinhStretch
from astropy.visualization.mpl_normalize import ImageNormalize
pl.matplotlib.rc_file('pubfiguresrc')

F = aplpy.FITSFigure(dpath("c18o_45to65_moment0.fits"), figure=pl.figure(1))
F.show_grayscale(invert=True, vmax=160)

F.recenter(290.92456, 14.51073, radius=0.04)
F.show_contour(dpath("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits"),
                colors=['r']*10, levels=np.logspace(-4,-2, 5), smooth=5) 

F.add_scalebar(length=((0.5 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.5 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)

F.save(fpath("diffuse/continuum_on_c18o.png"))
