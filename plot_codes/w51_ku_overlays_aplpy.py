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

distance = 5.41*u.kpc
e1e2 = coordinates.SkyCoord(290.93268,14.508363,unit=('deg','deg'), frame='icrs')

def set_tight_ticks(F):
    F.tick_labels.set_yformat('dd:mm:ss.ss')
    F.tick_labels.set_xformat('hh:mm:ss.ss')
    F.ticks.set_xspacing(0.001)
    F.ticks.set_yspacing(0.001)

#aplpy.make_rgb_cube( ('W51-CBAND-feathered.fits','W51-X-ABCD-S1.VTESS.VTC.DAVID-MEH.fits','W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'), 'W51_CXU_rgb' )

pl.close(1)
figure = pl.figure(1)
figure.clf()

# clean the header of junk axes
hdu = fits.open(dpath('W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'))
hdu[0].data = hdu[0].data.squeeze()
hdu[0].header = wcs.WCS(hdu[0].header).sub([wcs.WCSSUB_CELESTIAL]).to_header()

F = aplpy.FITSFigure(hdu,convention='calabretta',figure=figure)
#F = aplpy.FITSFigure(dpath+'W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.rot45.fits',convention='calabretta',figure=figure)
F.tick_labels.set_xformat('dd.dd')
F.tick_labels.set_yformat('dd.dd')
F.tick_labels.set_font(size=20)
F.axis_labels.set_font(size=20)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.011,invert=True)
#e1 = coordinates.ICRS(290.93263,14.50745,unit=('deg','deg'))
#F.recenter(e1.ra.value,e1.dec.value,width=1/60.,height=1/60.)
#F.recenter(290.92633,14.514769,radius=1.4/60.)
F.recenter(290.92345,14.511772,radius=1.5/60.)

F.add_scalebar(length=((0.5 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.5 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)

F.save(fpath('W51_Ku_grayscale.pdf'))

F.show_regions(rpath('HCHII_candidates.reg'))
#F.show_regions('/Users/adam/work/w51/cycle2_ALMA_frame2.reg')

F.save(fpath('W51_Ku_grayscale_HCHIIcandidates.pdf'))


F.scalebar.set_length(((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.1 pc')

F.recenter(290.91644,14.518939,radius=0.3/60.)
log.info("Reading briggs0_contsub image cube")
cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_briggs0_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
for velo in ProgressBar(np.arange(60,72,0.5)):
    #log.info("Velocity {0}".format(velo))
    c = pl.cm.jet_r((70-velo)/10.)
    #colors = [c[:3] + (x,) for x in (0.9,0.7,0.5,0.3,0.1)]
    #F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
    #               levels=[-1,-0.003,-0.002,-0.001],colors=colors,
    #               filled=False, layer='temporary')
    colors = [c[:3] + (x,) for x in (0.9,0.8,0.7,0.6,0.5)]
    colors = [c[:3] + (1,) for x in (0.9,0.8,0.7,0.6,0.5)]
    colors = ['r'] * 10
    hdu = cube[cube.closest_spectral_channel(velo*u.km/u.s),:,:].hdu
    if hdu.data.min() < -0.001:
        F.show_contour(hdu,
                       levels=[-1,-0.003,-0.002,-0.001],
                       colors=colors,
                       filled=False, layer='temporary')
        F.add_label(290.91254, 14.522828,
                    text="{0:0.1f} km s$^{{-1}}$".format(velo),
                    color='w', layer='label', size=20)
        F.recenter(290.91644,14.518939,radius=0.3/60.)
        F.save(fpath('contour_movie/IRS2_h2co22_on_cont22_v{0}.png'.format(velo)))

        F.add_label(290.93147, 14.508327,
                    text="{0:0.1f} km s$^{{-1}}$".format(velo),
                    color='w', layer='label', size=20)
        F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
        set_tight_ticks(F)
        F.save(fpath('contour_movie/e1e2_h2co22_on_cont22_briggs0_v{0}.png'.format(velo)))
        #log.info("Finished velo {0}".format(velo))
    else:
        log.info("No signal in {0}, continuing".format(velo))

    # not emission
    # F.show_regions(rpath('W51_22_emission_labels.reg'), layer='temporary2')
    # F.save(fpath('contour_movie/e1e2_h2co22_on_cont22_briggs0_v{0}_labeled.png'.format(velo)))
    # F.hide_layer('temporary2')

    F.hide_layer('temporary')

log.info("Reading natural_contsub image cube")
cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
for velo in ProgressBar(np.arange(60,72,0.5)):
    #log.info("Velocity {0}".format(velo))
    c = pl.cm.jet_r((70-velo)/10.)
    colors = [c[:3] + (x,) for x in (0.6,0.5,0.4,0.3,0.2,0.1)]
    colors = [c[:3] + (1,) for x in (0.6,0.5,0.4,0.3,0.2,0.1)]
    colors = ['r'] * 10
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[-1,-0.01,-0.005,-0.0025,-0.00125],colors=colors,
                   filled=False, layer='temporary')
    F.add_label(290.91254, 14.522828,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.recenter(290.91644,14.518939,radius=0.3/60.)
    F.save(fpath('contour_movie/IRS2_h2co22_on_cont22_natural_v{0}.png'.format(velo)))

    F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
    set_tight_ticks(F)
    F.add_label(290.93147, 14.508327,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.save(fpath('contour_movie/e1e2_h2co22_on_cont22_natural_v{0}.png'.format(velo)))

    # not emission
    # F.show_regions(rpath('W51_22_emission_labels.reg'), layer='temporary2')
    # F.save(fpath('contour_movie/e1e2_h2co22_on_cont22_natural_v{0}_labeled.png'.format(velo)))
    # F.hide_layer('temporary2')

    F.hide_layer('temporary')

cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
vr = 52,64
for velo in ProgressBar(np.arange(vr[0],vr[1]+0.5,0.5)):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.2,0.6,4)]
    colors = [c[:3] + (1,) for x in np.linspace(0.2,0.6,4)]
    colors = ['r'] * 10
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[0.002,0.004,0.006,1],colors=colors,
                   filled=False, layer='temporary')
    F.add_label(290.91254, 14.522828,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.recenter(290.91644,14.518939,radius=0.3/60.)
    F.save(fpath('contour_movie/IRS2_h2co22_emission_on_cont22_natural_v{0}.png'.format(velo)))

    F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
    set_tight_ticks(F)
    F.add_label(290.93147, 14.508327,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.save(fpath('contour_movie/e1e2_h2co22_emission_on_cont22_natural_v{0}.png'.format(velo)))
    F.show_regions(rpath('W51_22_emission_labels.reg'), layer='temporary2')
    F.save(fpath('contour_movie/e1e2_h2co22_emission_on_cont22_natural_v{0}_labeled.png'.format(velo)))
    F.remove_layer('temporary')
    F.remove_layer('temporary2')
    F.remove_layer('label')

#cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_briggs0_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
vr = 55,62
peak = cube.spectral_slab(vr[0]*u.km/u.s, vr[1]*u.km/u.s).max(axis=0)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.011,invert=True)
F.show_contour(peak.hdu,
               levels=[0.002, 0.004,0.006,0.008,1],
               colors=[(1,0,0,ii) for ii in np.linspace(0.5, 1, 5)],
               filled=False, layer='temporary')
F.show_regions(rpath('W51_22_emission_labels.reg'), layer='temporary2')
F.save(fpath('contour_movie/e1e2_h2co22_emission_on_cont22_natural_v{0}to{1}peak_labeled.png'.format(vr[0],vr[1])))
F.remove_layer('temporary')
F.remove_layer('temporary2')

F.recenter(290.91669,14.518151,radius=8./3600.)
F.scalebar.set_length(((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.05)
for velo in ProgressBar(np.arange(vr[0],vr[1]+0.5,0.5)):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.6,0.8,6)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[0.002,0.003,0.004,0.005,0.006,1],colors=colors,
                   linewidths=[1,2,3,4,5,6],
                   filled=False)
F.save(fpath('contour_movie/IRS2_h2co22_emission_on_cont22_natural_allvelos.png'))
F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
set_tight_ticks(F)
F.save(fpath('contour_movie/e1e2_h2co22_emission_on_cont22_natural_allvelos.png'))

for layer in F._layers.keys():
    if 'txt' not in layer: F.remove_layer(layer)

F.recenter(290.91669,14.518151,radius=8./3600.)
cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_briggs0_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio')
for velo in np.arange(vr[0],vr[1]+0.5,0.5):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.6,0.8,6)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[0.0015, 0.002,0.003,0.004,0.005,0.006,1],colors=colors,
                   linewidths=[1,2,3,4,5,6],
                   filled=False)
F.save(fpath('contour_movie/IRS2_h2co22_emission_on_cont22_briggs0_allvelos.png'))
F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
set_tight_ticks(F)
F.save(fpath('contour_movie/e1e2_h2co22_emission_on_cont22_briggs0_allvelos.png'))

for layer in F._layers.keys():
    if 'txt' not in layer: F.remove_layer(layer)





"""
#F.show_contour(dpath+'H2CO_22_Ku_D_tausummed_52to58.fits',levels=[2.0,3.0,4.0,10.0],colors=[(1,0,0,0.1),(1,0,0,0.2),(1,0,0,0.3)],filled=True,slices=[0])
#F.show_contour('../w51/H2CO_11_C_C_tausummed_42to61.fits',levels=[2.0,5.5],colors=[(1,0,0,0.4),(1,0,0,0.6)],filled=True,slices=[0])
#F.show_contour('../w51/H2CO_11_C_C_tausummed_63to67.fits',levels=[2.5,10.],colors=[(0,0,1,0.4),(0,0,1,0.6)],filled=True,slices=[0])
#F.show_contour('../w51_c_c/H2CO_11_C_C_tausummed_68to71.fits',levels=[0.6,5.0],colors=[(1,0,0,0.4),(1,0,0,0.6)],filled=True,slices=[0])
F.show_contour(dpath('H2CO_22_Ku_D_tausummed_52to58.fits'),levels=[0.5,2,15.0],colors=[(0.6,0.1,0.6,0.3),(0.6,0.1,0.6,0.5)],filled=True,slices=[0])
F.show_contour(dpath('H2CO_22_Ku_D_tausummed_63to67.fits'),levels=[0.9,  15.0],colors=[(0.0,0.5,1.0,0.3),(0.0,0.5,1.0,0.5)],filled=True,slices=[0])
F.show_contour(dpath('H2CO_22_Ku_D_tausummed_67to70.fits'),levels=[0.4,  15.0],colors=[(1.0,0.5,0.0,0.3),(1.0,0.5,0.0,0.5)],filled=True,slices=[0])
F.show_contour(dpath('H2CO_22_Ku_BD_small_tausummed_59to68.fits'),levels=[0.1,0.2,0.4,  15.0],colors=[(1.0,0.5,0.0,0.3),(1.0,0.5,0.0,0.5),(1.0,0.5,0.0,0.7)],filled=True,slices=[0])
#F.show_contour('../w51_c_c/H2CO_11_C_C_absorbsummed_42to61.fits',levels=[-0.2,-0.3,-0.4],colors=[(1,0,0,0.1),(1,0,0,0.2),(1,0,0,0.3)],filled=True,slices=[0])
#F.show_contour('../w51_c_c/H2CO_11_C_C_absorbsummed_63to67.fits',levels=[-0.2,-0.3,-0.4],colors=[(0,1,1,0.1),(0,1,1,0.2),(0,1,1,0.3)],filled=True,slices=[0])
#F.show_contour('../w51_c_c/H2CO_11_C_C_absorbsummed_68to71.fits',levels=[-0.2,-0.3,-0.4],colors=[(1,0,1,0.1),(1,0,1,0.2),(1,0,1,0.3)],filled=True,slices=[0])
F.show_contour(dpath+'v2.0_ds2_l050_13pca_map20.fits',levels=[12],colors=['b'], linewidths=[2], filled=False,convention='calabretta')

#F.show_contour(dpath+'SgrB2_nh3_3-3_maximum.fits',levels=[0.04],colors=['b'])
#F.show_contour(dpath+'SGRB2_1.3CM_fix_gal.fits',levels=[0.003],colors=['r'],smooth=3)
#F.show_rectangles(sgrb2m.l.value-0.01,sgrb2m.b.value-0.005,width=5/60.,height=5/60.,color='g')

#F.show_regions(dpath+'scalebars_8.5kpc_gal.reg')
#F.show_regions(dpath+'nh3_observed_region.reg')


beamd=(u.radian*((220*u.GHz).to(u.m,u.spectral()))/(12*u.m)).to(u.degree)
print beamd.to(u.arcsec)
beamd=beamd.value
#F.add_beam(major=beamd/2,minor=beamd/2,hatch='///',facecolor='none',color=(1,1,0,1))
#F.add_beam(major=0.5/3600.,minor=0.5/3600.,hatch='|||',facecolor='none',color=(0,1,1,1))


F.save(fpath('W51_Ku_withH2COcontours_noinset.pdf'))


F.hide_layer('contour_set_1')
F.hide_layer('contour_set_2')
F.save(fpath('W51_Ku_withH2COcontours_red.pdf'))
F.hide_layer('contour_set_3')
F.show_layer('contour_set_2')
F.save(fpath('W51_Ku_withH2COcontours_green.pdf'))
F.hide_layer('contour_set_2')
F.show_layer('contour_set_1')
F.save(fpath('W51_Ku_withH2COcontours_blue.pdf'))


(xl,yl),(xu,yu) = F._ax1.bbox._bbox.corners()[[0,3]]

inset = aplpy.FITSFigure(hdu,convention='calabretta',figure=figure,subplot=[xl,yu-0.25,0.2,0.25])
inset.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.031)
e1e2 = coordinates.SkyCoord(290.93268,14.508363,unit=('deg','deg'), frame='icrs')
inset.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
inset.tick_labels.hide()
inset.axis_labels.hide()
inset.show_contour(dpath('W51Ku_BD_spw19.bigish_uniform_contsub19.clean.image.integ_52to65.fits'), levels=[0.005,0.01], colors=['r','r'], smooth=3, slices=[0], linewidths=[2,2])
#inset.add_beam(major=0.5/3600., minor=0.5/3600., color='orange', linewidth=3, label='1"')
inset.add_scalebar(length=((0.1 * u.pc)/distance*u.radian).to(u.degree).value, color='orange', linewidth=3, label='0.1 pc')


F.save(fpath('W51_Ku_withH2COcontours.png'))
F.save(fpath('W51_Ku_withH2COcontours.pdf'))
"""
