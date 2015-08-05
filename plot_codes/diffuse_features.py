import numpy as np
import paths
from astropy.io import fits
import wcsaxes
from wcsaxes import WCS as WCSaxes
from astropy import wcs
import pylab as pl
from astropy.visualization import AsinhStretch,LinearStretch
from astropy.visualization.mpl_normalize import ImageNormalize
pl.matplotlib.rc_file('pubfiguresrc')

annotations = {'w51main_peak': [[(290.92755, 14.511689), [-3, 3], [-20, 20], 'r'],
                                [(290.92839, 14.510252), [-5, 0], [-30, 0], 'r'],
                                [(290.92643, 14.512817), [-10, 0], [-35, 0], 'b'],
                                [(290.92375, 14.511762), [-10, 0], [-40, 0], 'b'],
                                #[(290.92358, 14.510963), [10, -10], [30, -30], 'b'],
                                [(290.92064, 14.510635), [0, 0], [-52, 78], 'c'],
                               ],
               'irs2_C_low' : [[(290.91379, 14.522722), [-5, 5], [-20, 30], 'b'],
                               [(290.92188, 14.518806), [0, 0], [30, 0], 'r'],
                               [(290.92088, 14.516285), [0, 0], [20, 20], 'r'],
                               [(290.91549, 14.51498), [0, 0], [-5, 25], 'r'],
                              ],
              }
annotations['w51main_peak_diff'] = annotations['w51main_peak']
annotations['irs2_C_high'] = annotations['irs2_C_low'] + [[(290.91597, 14.51945), (0, -2), (0, -30), 'c']]
annotations['irs2_C_diff'] = annotations['irs2_C_low']
annotations['irs2_Ku_low'] = annotations['irs2_C_low']

def annotate(ax, xy, to_offset, from_offset, color, mywcs):

    ann = pl.Annotation(".", xy=mywcs.wcs_world2pix([xy],0)[0]+np.array(to_offset),
                        xytext=mywcs.wcs_world2pix([xy],0)[0]+np.array(from_offset),
                        arrowprops=dict(arrowstyle="->", edgecolor=color,
                                        color=color), color=color,
                        clip_on=False)
    ann.set_text("")
    ax.add_artist(ann)

for fn,pfx,coord_limits, (vmin,vmax), name, stretch in (
     ("W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits", 'Ku', [ (290.94225, 14.505832, ), (290.93794,  14.509374, ),], (-1e-1,3e-1), 'e11_bow', LinearStretch()),
     ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.9304, 14.5083), (290.9194, 14.5189)], (0.01, 9), 'w51main_peak', AsinhStretch()),
     ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.92462, 14.516962), (290.92352, 14.517891)], (2e-2, 2.0), 'e6', AsinhStretch()),
     ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.90023, 14.523703), (290.89873, 14.525156)], (2e-2, 0.3), 'd3', AsinhStretch()),
     ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.93729, 14.485868), (290.93594, 14.487230)], (-2e-2, 0.3), 'e7', AsinhStretch()),
     ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.93283, 14.506909), (290.93200, 14.507676)], (6e-2, 9), 'e1', AsinhStretch()),
     ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.92402, 14.513314), (290.90971, 14.525246)], (-6e-2, 5), 'irs2_C_low', AsinhStretch()),
     ("W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits", 'Ku', [ (290.92402, 14.513314), (290.90971, 14.525246)], (-6e-2, 5), 'irs2_Ku_low', AsinhStretch()),
     ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.92402, 14.513314), (290.90971, 14.525246)], (0, 13), 'irs2_C_high', AsinhStretch()),
     ("W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits", 'Ku', [ (290.92402, 14.513314), (290.90971, 14.525246)], (0, 52), 'irs2_Ku_high', AsinhStretch()),
     ("Cband_Epoch3sm-Epoch3.fits", 'C', [ (290.9304, 14.5083), (290.9194, 14.5189)], (-0.2, 1), 'w51main_peak_diff', LinearStretch()),
     ("Cband_Epoch3sm-Epoch3.fits", 'C', [ (290.92402, 14.513314), (290.90971, 14.525246)], (-0.2, 1.0), 'irs2_C_diff', LinearStretch()),
    ):

    hdu = fits.open(paths.dpath(fn))[0]
    mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
    wcsaxes = WCSaxes(mywcs.to_header())


    fig = pl.figure(1)
    fig.clf()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)

    im = ax.imshow(hdu.data.squeeze()*1e3, cmap=pl.cm.gray_r, origin='lower',
                   vmin=vmin, vmax=vmax,
                   norm=ImageNormalize(stretch=stretch))

    (x1,y1),(x2,y2) = mywcs.wcs_world2pix(coord_limits, 0)
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
    ra = ax.coords['ra']
    ra.set_major_formatter('hh:mm:ss.s')
    dec = ax.coords['dec']
    ra.set_axislabel("RA")
    dec.set_axislabel("Dec")
    cb = fig.colorbar(im)
    cb.set_label("mJy/beam")

    if name in annotations:
        for xy, to_offset, from_offset, color in annotations[name]:
            annotate(ax, xy, to_offset, from_offset, color, mywcs)

    pl.draw(); pl.show()
    pl.savefig(paths.fpath("diffuse/{0}.png".format(name)), bbox_inches='tight')

