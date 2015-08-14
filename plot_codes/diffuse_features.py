import numpy as np
from astropy import log
from astropy import units as u
import aplpy
import paths
from astropy.io import fits
import wcsaxes
from wcsaxes import WCS as WCSaxes
from astropy import wcs
import pylab as pl
from astropy.visualization import AsinhStretch,LinearStretch
from astropy.visualization.mpl_normalize import ImageNormalize
import common_constants
from common_constants import distance
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
     ("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits", 'C', [ (290.93024, 14.506951), (290.92676, 14.509907)], (-0.02, 1), 'peak_cluster_C_high', LinearStretch()),
     ("W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits", 'Ku', [ (290.92402, 14.513314), (290.90971, 14.525246)], (0, 52), 'irs2_Ku_high', AsinhStretch()),
     ("Cband_Epoch3sm-Epoch3.fits", 'C', [ (290.9304, 14.5083), (290.9194, 14.5189)], (-0.2, 1), 'w51main_peak_diff', LinearStretch()),
     ("Cband_Epoch3sm-Epoch3.fits", 'C', [ (290.92402, 14.513314), (290.90971, 14.525246)], (-0.2, 1.0), 'irs2_C_diff', LinearStretch()),
     ("Cband_Epoch3sm-Epoch3.fits", 'C', [ (290.93024, 14.506951), (290.92676, 14.509907)], (-0.02, 0.5), 'peak_cluster_C_diff', LinearStretch()),
    ):

    log.info("file {0} name {1}".format(fn, name))

    hdu = fits.open(paths.dpath(fn))[0]
    mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
    wcsaxes = WCSaxes(mywcs.to_header())


    fig = pl.figure(1, figsize=(10,10))
    fig.clf()

    F = aplpy.FITSFigure(hdu, subplot=[0.15, 0.1, 0.8, 0.8], figure=fig)

    F.show_grayscale(invert=True, vmin=vmin/1e3, vmax=vmax/1e3, stretch='linear' if
                     stretch == LinearStretch() else 'arcsinh')

    cx = (coord_limits[0][0] + coord_limits[1][0])/2.
    cy = (coord_limits[0][1] + coord_limits[1][1])/2.
    # radius = diamater / 2.
    coord_diffs = np.abs(np.array(coord_limits[0])-np.array(coord_limits[1]))
    radius = np.max(coord_diffs)/2.
    log.info("circle({0}, {1}, {2})".format(cx, cy, radius))
    F.recenter(cx, cy, radius=radius)

    if radius > 0.006:
        F.tick_labels.set_yformat('dd:mm:ss.s')
        F.tick_labels.set_xformat('hh:mm:ss.s')
        F.ticks.set_xspacing(0.005)
        F.ticks.set_yspacing(0.005)
    elif radius > 0.003:
        F.tick_labels.set_yformat('dd:mm:ss.ss')
        F.tick_labels.set_xformat('hh:mm:ss.ss')
        F.ticks.set_xspacing(0.001)
        F.ticks.set_yspacing(0.001)
    elif radius < 0.001:
        F.tick_labels.set_yformat('dd:mm:ss.sss')
        F.tick_labels.set_xformat('hh:mm:ss.sss')
        F.ticks.set_xspacing(0.0025)
        F.ticks.set_yspacing(0.0025)
    else:
        F.tick_labels.set_yformat('dd:mm:ss.ss')
        F.tick_labels.set_xformat('hh:mm:ss.ss')
        F.ticks.set_xspacing(0.001)
        F.ticks.set_yspacing(0.001)

    if radius > 0.003:
        sblength = 0.1*u.pc
    else:
        sblength = 0.05*u.pc
    F.add_scalebar(length=(sblength/distance*u.radian).to(u.degree).value)
    F.scalebar.set_label(sblength)
    F.scalebar.set_color('black')
    F.scalebar.set_linewidth(3)
    F.scalebar.set_font_size(20)
    pixscale = np.abs(F._wcs.pixel_scale_matrix[1,1])
    if name in annotations:
        for (x,y), to_offset, from_offset, color in annotations[name]:
            dx = (from_offset[0]-to_offset[0]) / pixscale
            dy = (from_offset[1]-to_offset[1]) / pixscale
            F.show_arrows(x, y, dx, dy, color=color)

    F.show_colorbar()
    F.colorbar.set_axis_label_text("mJy/beam")
    F.save(paths.fpath("diffuse/{0}_aplpy.png".format(name)))

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

    sblength_deg = (sblength/distance).to(u.degree, u.dimensionless_angles())
    sblength_pix = sblength_deg.value * pixscale
    ax.plot([x1 + np.abs(x2-x1)*0.05,
             x1 + np.abs(x2-x1)*0.05 + sblength_pix],
            [y1 + np.abs(y2-y1)*0.05]*2,
            linewidth=3,
            color='black')
    ax.text(np.mean([x1 + np.abs(x2-x1)*0.05,
                     x1 + np.abs(x2-x1)*0.05 + sblength_pix]),
            y1 + np.abs(y2-y1)*0.05,
            s=str(sblength))

    if name in annotations:
        for xy, to_offset, from_offset, color in annotations[name]:
            annotate(ax, xy, to_offset, from_offset, color, mywcs)
    pl.draw(); pl.show()
    pl.savefig(paths.fpath("diffuse/{0}.png".format(name)), bbox_inches='tight')
