from astropy import coordinates
from astropy import units as u
from astropy.io import fits
from astropy import wcs
import aplpy
import pylab
import paths

dpath = '/Volumes/128gbdisk/w51/'
dpath = '/Users/adam/work/w51/paper_w51_evla/data/'

#aplpy.make_rgb_cube( ('W51-CBAND-feathered.fits','W51-X-ABCD-S1.VTESS.VTC.DAVID-MEH.fits','W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'), 'W51_CXU_rgb' )

for fn,(vmin,vmax),name in (("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits",(-0.6e-3,13e-3), 'W51_C'),
                            ('W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits',(-2.0e-3,45e-3), 'W51_Ku')):
    figure = pylab.figure(1)
    figure.clf()

    # clean the header of junk axes
    hdu = fits.open(dpath+fn)
    for ii in [3,4]:
        for kk in ['CRVAL','CTYPE','CDELT','CRPIX','CUNIT','NAXIS']:
            k = kk+str(ii)
            if k in hdu[0].header:
                del hdu[0].header[k]
        for jj in [1,2,3,4]:
            k = "PC{0:02d}_{1:02d}".format(ii,jj)
            if k in hdu[0].header:
                del hdu[0].header[k]
            k = "PC{1:02d}_{0:02d}".format(ii,jj)
            if k in hdu[0].header:
                del hdu[0].header[k]

    hdu[0].header['NAXIS'] = 2
    hdu[0].data = hdu[0].data.squeeze()
    w = wcs.WCS(hdu[0].header)
    assert w.naxis == 2

    F = aplpy.FITSFigure(hdu,convention='calabretta',figure=figure)
    #F = aplpy.FITSFigure(dpath+'W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.rot45.fits',convention='calabretta',figure=figure)
    #F.tick_labels.set_xformat('dd.dd')
    #F.tick_labels.set_yformat('dd.dd')
    F.tick_labels.set_font(size=20)
    F.axis_labels.set_font(size=20)
    F.show_grayscale(stretch='arcsinh',vmin=vmin,vmax=vmax, invert=True)
    #e1 = coordinates.ICRS(290.93263,14.50745,unit=('deg','deg'))
    #F.recenter(e1.ra.value,e1.dec.value,width=1/60.,height=1/60.)
    #F.recenter(290.92633,14.514769,radius=1.4/60.)
    F.recenter(290.92345,14.511772,radius=1.1/60.)
    F.add_scalebar(length=((0.5 * u.pc)/(5.4*u.kpc)*u.radian).to(u.degree).value)
    F.scalebar.set_label('0.5 pc')
    F.scalebar.set_color('orange')
    F.scalebar.set_linewidth(3)
    F.scalebar.set_font_size(20)
    F.save(paths.fpath('diffuse/{0}_overview.pdf'.format(name)), dpi=300)

    F.show_regions(paths.rpath("diffuse_hii_regions.reg"))
    F.save(paths.fpath('diffuse/{0}_overview_diffusehiiregionlabels.pdf'.format(name)), dpi=300)

    F.remove_layer('region_set_1')
    F.show_regions(paths.rpath("pointsource_centroids.reg"))
    F.save(paths.fpath('diffuse/{0}_overview_pointsourcelabels.pdf'.format(name)), dpi=300)
