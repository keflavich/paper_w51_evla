import warnings
import numpy as np
from astropy import units as u
from spectral_cube import SpectralCube
from astropy.io import fits
import paths
from paths import dpath,rpath,fpath,pcpath
import glob
import pyspeckit
import os
import pyregion
import time
import fnames
import matplotlib
from astropy.table import Table,Column
matplotlib.rc_file(pcpath('pubfiguresrc'))
import radio_beam

import pylab as pl
pl.ioff()

cube_name_titles = {'11_natural': 'H$_2$CO 1-1 natural',
                    '22_natural': 'H$_2$CO 2-2 natural',
                    '22_briggs0': 'H$_2$CO 2-2 Robust 0',
                    '11_uniform': 'H$_2$CO 1-1 uniform',}


for cube_name, cube_fn in fnames.cube_names.items():

    cube = SpectralCube.read(dpath(cube_fn)).with_spectral_unit(u.km/u.s,
                                                                velocity_convention='radio')
    errspec = cube.std(axis=(1,2))
    pix_area = np.abs(np.product(np.diag(cube.wcs.celestial.pixel_scale_matrix)))
    beam = radio_beam.Beam.from_fits_header(cube.header)
    ppbeam = beam.sr.to(u.deg**2).value/pix_area

    JyToK = u.Jy.to(u.K, equivalencies=u.brightness_temperature(beam,
                                                                cube.wcs.wcs.restfrq*u.Hz))

    for regname, outdir in [('W51_22_emission.reg', 'emission'),
                            ('W51_e_apertures.reg', 'hiiregionh2co')]:

        outpath = dpath('spectra/{0}'.format(outdir))
        figpath = fpath('spectra/{0}'.format(outdir))

        if not os.path.exists(outpath):
            os.makedirs(outpath)
        if not os.path.exists(figpath):
            os.makedirs(figpath)


        regions = pyregion.open(rpath(regname))


        prefix = os.path.basename(cube_fn)

        t0 = time.time()
        
        for R in regions:
            subcube = cube.subcube_from_ds9region(pyregion.ShapeList([R]))
            pcube = pyspeckit.Cube(cube=subcube)

            if 'text' not in R.attr[1]:
                print "Skipping region {0}".format(R)
                continue
            elif R.name not in ('circle','ellipse'):
                print "Skipping region {0} because it has the wrong shape".format(R)
                continue


            name = R.attr[1]['text']
            print "dt=%g" % (time.time()-t0), name,R
            sp = pcube.get_apspec(R.coord_list,coordsys='celestial',
                                  wunit='degree', method='mean')
            if not np.any(sp.data):
                warnings.warn("No data found for {0} in {1}".format(name, cube_fn))
                continue
            sp.specname = name
            #sp.data /= ppbeam
            sp.error = (sp[:sp.xarr.x_to_pix(45)].stats()['std']*errspec/errspec.min()).value
            #sp.error = errspec/ppbeam**0.5
            if 'APRADIUS' in sp.header:
                sp.header['APRADPIX'] = sp.header['APRADIUS'] / np.abs(sp.header['OLDCDEL1']*sp.header['CDELT2'])**0.5
                sp.header['APAREA'] = sp.header['APRADIUS']**2 * np.pi
            else:
                sp.header['APRADPIX'] = (sp.header['APMAJ'] *
                                         sp.header['APMIN'] /
                                         np.abs(sp.header['OLDCDEL1'] *
                                                sp.header['CDELT2']))**0.5
                sp.header['APAREA'] = sp.header['APREFF']**2 * np.pi
            sp.header['APAREAPX'] = sp.header['APRADPIX']**2 * np.pi

            sp.header['FILENAME'] = cube_fn
            sp.header['JYTOK'] = JyToK
            sp.write(os.path.join(outpath, '%s_%s.fits' % (prefix,name)))

            pl.figure(1).clf()
            sp.plotter(errstyle='fill', figure=pl.figure(1))
            sp.plotter.axis.set_title("{0}: {1}".format(cube_name_titles[cube_name], name))
            sp.plotter.figure.savefig(os.path.join(figpath, '%s_%s.png' %
                                                   (prefix,name)))

            # Jy -> mJy
            sp.data *= 1e3
            sp.error *= 1e3
            #sp.unit = '$T_B$ (K)'
            sp.unit = 'mJy/beam'
            sp.plotter(errstyle='fill')
            ax2 = sp.plotter.axis.twinx()
            ax2.set_ylim(*(np.array(sp.plotter.axis.get_ylim()) * JyToK/1e3))
            ax2.set_ylabel("$T_B$ (K)")

            sp.plotter.figure.savefig(os.path.join(figpath, '%s_%s_K.png' %
                                                   (prefix,name)))

pl.ion()
