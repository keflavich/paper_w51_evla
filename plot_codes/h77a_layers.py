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
import warnings


regions = pyregion.open(paths.rpath('RRL_subframes.reg'))
cube = SpectralCube.read(paths.dpath('H77a_BDarray_speccube_briggs0_contsub_cvel_big.fits'))

for reg in regions:
    name = reg.attr[1]['text']
    reg = pyregion.ShapeList([reg])

    scube = cube.subcube_from_ds9region(reg)

    layer1 = scube.spectral_slab(22*u.km/u.s, 44*u.km/u.s).moment0()
    layer2 = scube.spectral_slab(44*u.km/u.s, 66*u.km/u.s).moment0()
    layer3 = scube.spectral_slab(66*u.km/u.s, 88*u.km/u.s).moment0()
    total = scube.spectral_slab(26*u.km/u.s, 86*u.km/u.s).moment0()

    layer1.hdu.writeto(paths.dpath('rrl_moments/{0}_layer1.fits'.format(name)), clobber=True)
    layer2.hdu.writeto(paths.dpath('rrl_moments/{0}_layer2.fits'.format(name)), clobber=True)
    layer3.hdu.writeto(paths.dpath('rrl_moments/{0}_layer3.fits'.format(name)), clobber=True)
    total.hdu.writeto(paths.dpath('rrl_moments/{0}_total.fits'.format(name)), clobber=True)

    F = aplpy.FITSFigure(total.hdu)
    F.show_grayscale(invert=True)

    F.show_contour(layer1.hdu, smooth=3, levels=np.linspace(0.003, 0.12, 6), colors=['b']*10)
    F.show_contour(layer2.hdu, smooth=3, levels=np.linspace(0.003, 0.12, 6), colors=['g']*10)
    F.show_contour(layer3.hdu, smooth=3, levels=np.linspace(0.003, 0.12, 6), colors=['r']*10)
    F.save(paths.fpath('rrls/h77a_3layers_{0}.png'.format(name)))

    # Begin channel map code here
    pl.figure(3, figsize=(16,16)).clf()
    Nrows=Ncols=4
    fig, ax = pl.subplots(Nrows, Ncols,
                          sharex=True,
                          sharey=True, num=3)
    
    # integrate over velocities to make channel maps of a set width
    vstart = 15 #km/s
    vend = 95
    vstep = 5
    if int(vend-vstart)/vend > Nrows*Ncols:
        warnings.warn("More slices than rows*columns!")
    if int(vend-vstart)/vend < Nrows*Ncols:
        warnings.warn("Fewer slices than rows*columns")

    layers = [scube.spectral_slab(vv*u.km/u.s,
                                  (vv+5)*u.km/u.s).moment0() for ii,vv in enumerate(np.arange(vstart,vend,vstep))]
    # Determine the maximum value to display
    mx = np.max([np.nanmax(x).value for x in layers])


    for ii,vv in enumerate(np.arange(vstart, vend, vstep, dtype='int')):
        v1 = vv
        v2 = (vv+5)
        #pl.subplot(4, 4, ii+1)
        layer = layers[ii]
        im = ax[ii / 4, ii % 4].imshow(layer.value,
                                       norm=ImageNormalize(vmin=-0.001,
                                                           vmax=mx,
                                                           stretch=AsinhStretch()),
                                       cmap=pl.cm.gray_r)
        ax[ii / 4, ii % 4].annotate("${0:d} < v < {1:d}$".format(v1,v2), (0.1, 0.8),
                                    xycoords='axes fraction', color='k', fontsize='large')

    #fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([1.00, 0.05, 0.04, 0.9])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label("mJy km s$^{-1}$")

    pl.subplots_adjust(hspace=0,
                       wspace=0)

    for i in range(Nrows):
        for j in range(Ncols):
            if i == 0:
                ax[i,j].xaxis.set_ticks_position('top')
                pl.setp(ax[i,j].get_xticklabels(), visible=False)
                ax[i,j].xaxis.set_ticklabels([])
            elif i == Nrows-1:
                ax[i,j].xaxis.set_ticks_position('bottom')
                pl.setp(ax[i,j].get_xticklabels(), visible=True)
            else:
                ax[i,j].xaxis.set_ticks_position('none')
                pl.setp(ax[i,j].get_xticklabels(), visible=False)
                ax[i,j].xaxis.set_ticklabels([])

            if j == 0:
                ax[i,j].yaxis.set_ticks_position('left')
            elif j == Ncols-1:
                ax[i,j].yaxis.set_ticks_position('right')
                pl.setp(ax[i,j].get_yticklabels(), visible=False)
                ax[i,j].yaxis.set_ticklabels([])
            else:
                ax[i,j].yaxis.set_ticks_position('none')
                pl.setp(ax[i,j].get_yticklabels(), visible=False)
                ax[i,j].yaxis.set_ticklabels([])

    pl.subplots_adjust(hspace=0,
                       wspace=0)
    pl.savefig(paths.fpath('rrls/h77a_channelmaps_{0}.png'.format(name)), bbox_inches='tight')
    pl.savefig(paths.fpath('rrls/h77a_channelmaps_{0}_hires.png'.format(name)), bbox_inches='tight', dpi=300)
