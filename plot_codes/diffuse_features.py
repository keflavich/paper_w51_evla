import paths
from astropy.io import fits
import wcsaxes
from wcsaxes import WCS as WCSaxes
from astropy import wcs
import pylab as pl
from astropy.visualization import AsinhStretch,LinearStretch
from astropy.visualization.mpl_normalize import ImageNormalize
pl.matplotlib.rc_file('pubfiguresrc')

for fn,pfx,coord_limits, (vmin,vmax), name, stretch in (
    ("W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits", 'Ku', [ (290.94225, 14.505832, ), (290.93794,  14.509374, ),], (-2e-1,3e-1), 'e11_bow', AsinhStretch()),
    ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.92977, 14.508804), (290.92003, 14.515556)], (0.01, 9), 'w51main_peak', AsinhStretch()),
    ("Cband_Epoch3-Epoch1.fits", 'C', [ (290.92977, 14.508804), (290.92003, 14.515556)], (-0.6, 1), 'w51main_peak_diff', LinearStretch())):

    hdu = fits.open(paths.dpath(fn))[0]
    mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
    wcsaxes = WCSaxes(mywcs.to_header())


    fig = pl.figure(1)
    fig.clf()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)

    im = ax.imshow(hdu.data.squeeze()*1e3, cmap=pl.cm.gray_r, origin='lower',
                   vmin=vmin, vmax=vmax,
                   norm=ImageNormalize(stretch=stretch))

    (x1,y1),(x2,y2) = mywcs.wcs_world2pix( coord_limits, 0)
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    ra.set_axislabel("RA")
    dec.set_axislabel("Dec")
    cb = fig.colorbar(im)
    cb.set_label("mJy/beam")

    pl.draw(); pl.show()
    pl.savefig(paths.fpath("diffuse/{0}.png".format(name)), bbox_inches='tight')

