"""
TODO: Re-do this on a larger (in velocity) H77a cube to see if
(1) there is a HeII counterpart (probably not; too faint)
(2) there is a redshifted counterpart

"""
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
import paths

from pvextractor import extract_pv_slice
from pvextractor.geometry import Path

jet_endpoints = coordinates.SkyCoord([290.9186, 290.91424]*u.deg, [14.518876, 14.517977]*u.deg, frame='fk5')

xy = Path(jet_endpoints, width=2.8*u.arcsec)
for ii,fn in enumerate(('w51.neii.square.fits',
                        'w51.siv.square.fits',
                        'W51north_H77_Outflow_cutout.fits',
                        'H77a_cutout_1.fits',
                        'H77a_cutout_2.fits',
                       )):
    cube = SpectralCube.read(paths.dpath(fn)).with_spectral_unit(u.km/u.s,
                                                                 velocity_convention='radio',
                                                                 rest_value=14128610000.0*u.Hz).spectral_slab(-200*u.km/u.s, 200*u.km/u.s)
    pv = extract_pv_slice(cube.hdu, xy)
    fig = pl.figure(ii)
    fig.clf()
    FF = aplpy.FITSFigure(pv, figure=fig)
    FF.show_grayscale(aspect='auto')#, stretch='arcsinh')
    FF.recenter(0.0021, 0, width=0.0042, height=200e3)
    FF.save(paths.fpath('jetpv/{0}.png'.format(fn[:-5])))

