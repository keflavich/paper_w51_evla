import numpy as np
import pylab as pl
from paths import dpath,fpath,rpath
from spectral_cube import SpectralCube
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
import aplpy
from fnames import fnames, cube_names, continua

distance = 5.1*u.kpc
e1e2 = coordinates.ICRS(290.93268,14.508363,unit=('deg','deg'))

#aplpy.make_rgb_cube( ('W51-CBAND-feathered.fits','W51-X-ABCD-S1.VTESS.VTC.DAVID-MEH.fits','W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'), 'W51_CXU_rgb' )

pl.close(1)
figure = pl.figure(1)
figure.clf()

# clean the header of junk axes
hdu = fits.open(dpath(continua['11_cont_both_2048']))
for ii in [3,4]:
    for kk in ['CRVAL','CTYPE','CDELT','CRPIX','CUNIT','NAXIS']:
        k = kk+str(ii)
        if k in hdu[0].header:
            del hdu[0].header[k]
hdu[0].header['NAXIS'] = 2
hdu[0].data = hdu[0].data.squeeze()

F = aplpy.FITSFigure(hdu,convention='calabretta',figure=figure)
F.tick_labels.set_xformat('dd.dd')
F.tick_labels.set_yformat('dd.dd')
F.tick_labels.set_font(size=20)
F.axis_labels.set_font(size=20)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.011)
#e1 = coordinates.ICRS(290.93263,14.50745,unit=('deg','deg'))
#F.recenter(e1.ra.value,e1.dec.value,width=1/60.,height=1/60.)
#F.recenter(290.92633,14.514769,radius=1.4/60.)
F.recenter(290.92345,14.511772,radius=1.5/60.)

F.add_scalebar(length=((0.5 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.5 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)

F.save(fpath('W51_C_grayscale.pdf'))

F.show_regions(rpath('HCHII_candidates.reg'))
#F.show_regions('/Users/adam/work/w51/cycle2_ALMA_frame2.reg')

F.save(fpath('W51_C_grayscale_HCHIIcandidates.pdf'))


F.scalebar.set_length(((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.1 pc')

F.recenter(290.91644,14.518939,radius=0.3/60.)
cube = SpectralCube.read(dpath(cube_names['11_uniform'])).with_spectral_unit(u.km/u.s, velocity_convention='radio')
for velo in np.arange(60,72,0.5):
    c = pl.cm.jet_r((70-velo)/10.)
    colors = [c[:3] + (x,) for x in (0.9,0.7,0.5,0.3,0.1)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[-1,-0.003,-0.002,-0.001],colors=colors,
                   filled=True, layer='temporary')
    F.add_label(290.91254, 14.522828,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.recenter(290.91644,14.518939,radius=0.3/60.)
    F.save(fpath('contour_movie/IRS2_h2co11_on_cont11_v{0}.png'.format(velo)))

    F.add_label(290.93147, 14.508327,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
    F.save(fpath('contour_movie/e1e2_h2co11_on_cont11_briggs0_v{0}.png'.format(velo)))
    F.hide_layer('temporary')

cube = SpectralCube.read(dpath(cube_names['11_natural'])).with_spectral_unit(u.km/u.s, velocity_convention='radio')
for velo in np.arange(60,72,0.5):
    c = pl.cm.jet_r((70-velo)/10.)
    colors = [c[:3] + (x,) for x in (0.6,0.5,0.4,0.3,0.2,0.1)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[-1,-0.01,-0.005,-0.0025,-0.00125],colors=colors,
                   filled=True, layer='temporary')
    F.add_label(290.91254, 14.522828,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.recenter(290.91644,14.518939,radius=0.3/60.)
    F.save(fpath('contour_movie/IRS2_h2co11_on_cont11_natural_v{0}.png'.format(velo)))

    F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
    F.add_label(290.93147, 14.508327,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.save(fpath('contour_movie/e1e2_h2co11_on_cont11_natural_v{0}.png'.format(velo)))
    F.hide_layer('temporary')

cube = SpectralCube.read(dpath(cube_names['11_natural'])).with_spectral_unit(u.km/u.s, velocity_convention='radio')
vr = 52,64
for velo in np.arange(vr[0],vr[1]+0.5,0.5):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.2,0.6,4)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[0.002,0.004,0.006,1],colors=colors,
                   filled=True, layer='temporary')
    F.add_label(290.91254, 14.522828,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.recenter(290.91644,14.518939,radius=0.3/60.)
    F.save(fpath('contour_movie/IRS2_h2co11_emission_on_cont11_natural_v{0}.png'.format(velo)))

    F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
    F.add_label(290.93147, 14.508327,
                text="{0:0.1f} km s$^{{-1}}$".format(velo),
                color='w', layer='label', size=20)
    F.save(fpath('contour_movie/e1e2_h2co11_emission_on_cont11_natural_v{0}.png'.format(velo)))
    F.remove_layer('temporary')
    F.remove_layer('label')

F.recenter(290.91669,14.518151,radius=8./3600.)
F.scalebar.set_length(((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.05)
for velo in np.arange(vr[0],vr[1]+0.5,0.5):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.6,0.8,6)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[0.002,0.003,0.004,0.005,0.006,1],colors=colors,
                   linewidths=[1,2,3,4,5,6],
                   filled=False)
F.save(fpath('contour_movie/IRS2_h2co11_emission_on_cont11_natural_allvelos.png'))
F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
F.save(fpath('contour_movie/e1e2_h2co11_emission_on_cont11_natural_allvelos.png'))

for layer in F._layers.keys():
    if 'txt' not in layer: F.remove_layer(layer)

F.recenter(290.91669,14.518151,radius=8./3600.)
cube = SpectralCube.read(dpath(cube_names['11_uniform'])).with_spectral_unit(u.km/u.s, velocity_convention='radio')
for velo in np.arange(vr[0],vr[1]+0.5,0.5):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.6,0.8,6)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[0.0015, 0.002,0.003,0.004,0.005,0.006,1],colors=colors,
                   linewidths=[1,2,3,4,5,6],
                   filled=False)
F.save(fpath('contour_movie/IRS2_h2co11_emission_on_cont11_briggs0_allvelos.png'))
F.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
F.save(fpath('contour_movie/e1e2_h2co11_emission_on_cont11_briggs0_allvelos.png'))

for layer in F._layers.keys():
    if 'txt' not in layer: F.remove_layer(layer)
