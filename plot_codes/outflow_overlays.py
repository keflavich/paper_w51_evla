import numpy as np
import pylab as pl
from paths import dpath,fpath,rpath
from spectral_cube import SpectralCube
from astropy import coordinates
from astropy import units as u
from astropy.io import fits
import aplpy
import paths
from common_constants import distance
import FITS_tools

e1e2 = coordinates.ICRS(290.93268,14.508363,unit=('deg','deg'))

#aplpy.make_rgb_cube( ('W51-CBAND-feathered.fits','W51-X-ABCD-S1.VTESS.VTC.DAVID-MEH.fits','W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'), 'W51_CXU_rgb' )

pl.close(1)
figure = pl.figure(1)
figure.clf()

# clean the header of junk axes
cont22hdu = fits.open(dpath('W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'))
for ii in [3,4]:
    for kk in ['CRVAL','CTYPE','CDELT','CRPIX','CUNIT','NAXIS']:
        k = kk+str(ii)
        if k in cont22hdu[0].header:
            del cont22hdu[0].header[k]
cont22hdu[0].header['NAXIS'] = 2
cont22hdu[0].data = cont22hdu[0].data.squeeze()

F = aplpy.FITSFigure(cont22hdu,convention='calabretta',figure=figure)
F.set_auto_refresh(False)
#F = aplpy.FITSFigure(dpath+'W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.rot45.fits',convention='calabretta',figure=figure)
F.tick_labels.set_xformat('dd.ddd')
F.tick_labels.set_yformat('dd.ddd')
F.tick_labels.set_font(size=20)
F.axis_labels.set_font(size=20)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.011)

F.add_scalebar(length=((0.5 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.5 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)


siv = SpectralCube.read(paths.dpath('w51.siv.fits'))
neii = SpectralCube.read(paths.dpath('w51.neii.fits'))
h77a = SpectralCube.read(paths.dpath('W51north_H77_Outflow_cutout.fits'))

siv_outflow = siv.spectral_slab(-28*u.km/u.s, -70*u.km/u.s).sum(axis=0)
neii_outflow = neii.spectral_slab(-28*u.km/u.s, -70*u.km/u.s).sum(axis=0)
h77a_outflow = h77a.spectral_slab(-16*u.km/u.s, -60*u.km/u.s).sum(axis=0)

c = (1,0.4,0)
sivcolors = [c[:3] + (x,) for x in (0.7,0.6,0.5,0.4,0.3,0.2,0.1)[::-1]]

F.recenter(290.91669,14.518151,radius=8./3600.)
F.scalebar.set_length(((0.1 * u.pc)/distance*u.radian).to(u.degree).value)
F.scalebar.set_label('0.1 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.05)

F.save(fpath('IRS2_W51_Ku_grayscale.pdf'), dpi=150)
F.save(fpath('IRS2_W51_Ku_grayscale.png'), dpi=150)

F.show_contour(siv_outflow.hdu, levels=[20,30,40,50,60,70,80],colors=sivcolors,
               filled=True, layer='temporary')
F.save(fpath('irs2outflow/IRS2_siv_on_cont22.png'), dpi=150)

c = (0,0.4,1)
neiicolors = [c[:3] + (x,) for x in (0.93,0.92,0.91,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)[::-1]]
F.show_contour(neii_outflow.hdu, levels=np.arange(30,190,10), colors=neiicolors,
               filled=True, layer='temporary')
F.save(fpath('irs2outflow/IRS2_neii_on_cont22.png'), dpi=150)

c = (0,0.9,0.1)
h77acolors = [c[:3] + (x,) for x in (0.1,0.2,0.3,0.4,0.5,0.6,0.7)]
h77alevels = np.arange(0.01,0.05,0.005)
F.show_contour(h77a_outflow.hdu, levels=h77alevels, colors=h77acolors,
               filled=True, layer='temporary')
F.save(fpath('irs2outflow/IRS2_h77a_on_cont22.png'), dpi=150)

F.remove_layer('temporary')

F.show_contour(neii_outflow.hdu, levels=np.arange(30,190,10), colors=neiicolors,
               filled=True, layer='neii')
F.show_contour(siv_outflow.hdu, levels=[20,30,40,50,60,70,80],colors=sivcolors,
               filled=True, layer='siv')

cutout_coords = {'xlo':290.91974*u.deg, 'xhi':290.90926*u.deg, 'ylo':14.51492*u.deg, 'yhi':14.523815*u.deg}
cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio').subcube(**cutout_coords)


F.show_contour(neii_outflow.hdu, levels=np.arange(30,190,10), colors=neiicolors,
               filled=True, layer='neii')
F.show_contour(siv_outflow.hdu, levels=[20,30,40,50,60,70,80],colors=sivcolors,
               filled=True, layer='siv')

cutout_coords = {'xlo':290.91974*u.deg, 'xhi':290.90926*u.deg, 'ylo':14.51492*u.deg, 'yhi':14.523815*u.deg}
cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_natural_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio').subcube(**cutout_coords)

vr = [56,60]
for velo in np.arange(vr[0],vr[1]+0.5,0.5):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.8,0.9,6)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[0.003,0.005,0.007,1],colors=colors,
                   filled=False,
                   linewidths=[1,2,3,4],
                   layer='h2co_{0}'.format(velo))

F.save(fpath('irs2outflow/IRS2_core_and_siv_and_neii_on_cont22.png'), dpi=150)
F.show_contour(h77a_outflow.hdu, levels=h77alevels,colors=h77acolors,
               filled=True, layer='h77a')
F.save(fpath('irs2outflow/IRS2_core_and_siv_and_neii_and_h77a_on_cont22.png'), dpi=150)

for layer in F._layers.keys():
    F.remove_layer(layer)


vr = [56,60]
for velo in np.arange(vr[0],vr[1]+0.5,0.5):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.8,0.9,6)]
    F.show_contour(cube[cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[0.003,0.005,0.007,1],colors=colors,
                   filled=False,
                   linewidths=[1,2,3,4],
                   layer='h2co_{0}'.format(velo))

F.save(fpath('irs2outflow/IRS2_core_on_cont22.png'), dpi=150)

for layer in F._layers.keys():
    F.remove_layer(layer)



briggs_h2co22cube = SpectralCube.read(dpath('W51Ku_BD_h2co_v30to90_briggs0_contsub.image.fits')).with_spectral_unit(u.km/u.s, velocity_convention='radio').subcube(**cutout_coords)
vr = [60,71]
h2co22_integrated = briggs_h2co22cube.spectral_slab(vr[0]*u.km/u.s, vr[1]*u.km/u.s).sum(axis=0)

cont_regrid = FITS_tools.hcongrid.hcongrid_hdu(cont22hdu[0], h2co22_integrated.header)
# out = (exp(-tau))*in
npix = briggs_h2co22cube.spectral_slab(vr[0]*u.km/u.s, vr[1]*u.km/u.s).shape[0]
tau = -np.log((cont_regrid.data*npix+h2co22_integrated.value) / (cont_regrid.data*npix))
tau[h2co22_integrated.value > -0.005] = np.nan
tau_hdu = fits.PrimaryHDU(data=tau, header=cont_regrid.header)

levels = np.logspace(np.log10(0.005),np.log10(0.2),6)
F.show_contour(tau_hdu,
               levels=levels,
               colors=[(0,0.3,0.9,x) for x in np.linspace(0.15,0.7,len(levels))],
               filled=True,
               #linewidths=#np.linspace(3,1,len(levels)),
               layer='h2co_integ')

F.save(fpath('irs2outflow/IRS2_opticaldepth_absorptioncontours_on_cont22_integrated.png'), dpi=150)

levels = [-0.04,-0.03,-0.02,-0.01]
F.show_contour(h2co22_integrated.hdu,
               levels=levels,
               colors=[(0,0.3,0.9,0.9)]*10,
               filled=False,
               linewidths=np.linspace(3,1,len(levels)),
               layer='h2co_integ')

F.save(fpath('irs2outflow/IRS2_absorptioncontours_on_cont22_integrated.png'), dpi=150)

F.remove_layer('h2co_integ')

for velo in np.arange(vr[0],vr[1]+0.5,0.5):
    c = pl.cm.jet_r((vr[1]-velo)/(vr[1]-vr[0]))
    colors = [c[:3] + (x,) for x in np.linspace(0.4,0.9,6)][::-1]
    F.show_contour(briggs_h2co22cube[briggs_h2co22cube.closest_spectral_channel(velo*u.km/u.s)].hdu,
                   levels=[-0.005,-0.004,-0.003,-0.002],colors=colors,
                   filled=False,
                   linewidths=[1,2,3,4,5][::-1],
                   layer='h2co_{0}'.format(velo))

F.save(fpath('irs2outflow/IRS2_absorptioncontours_on_cont22_velocities.png'), dpi=150)
