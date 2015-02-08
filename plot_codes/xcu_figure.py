from astropy import coordinates
from astropy import units as u
from astropy.io import fits
import aplpy
import pylab

dpath = '/Volumes/128gbdisk/w51/'

#aplpy.make_rgb_cube( ('W51-CBAND-feathered.fits','W51-X-ABCD-S1.VTESS.VTC.DAVID-MEH.fits','W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits'), 'W51_CXU_rgb' )

figure = pylab.figure(1)
figure.clf()

# clean the header of junk axes
hdu = fits.open(dpath+'W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits')
for ii in [3,4]:
    for kk in ['CRVAL','CTYPE','CDELT','CRPIX','CUNIT','NAXIS']:
        k = kk+str(ii)
        if k in hdu[0].header:
            del hdu[0].header[k]
hdu[0].header['NAXIS'] = 2
hdu[0].data = hdu[0].data.squeeze()

F = aplpy.FITSFigure(hdu,convention='calabretta',figure=figure)
#F = aplpy.FITSFigure(dpath+'W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.rot45.fits',convention='calabretta',figure=figure)
F.tick_labels.set_xformat('dd.dd')
F.tick_labels.set_yformat('dd.dd')
F.tick_labels.set_font(size=20)
F.axis_labels.set_font(size=20)
F.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.011)
#e1 = coordinates.ICRS(290.93263,14.50745,unit=('deg','deg'))
#F.recenter(e1.ra.value,e1.dec.value,width=1/60.,height=1/60.)
#F.recenter(290.92633,14.514769,radius=1.4/60.)
F.recenter(290.92345,14.511772,radius=1.1/60.)

#F.show_contour(dpath+'H2CO_22_Ku_D_tausummed_52to58.fits',levels=[2.0,3.0,4.0,10.0],colors=[(1,0,0,0.1),(1,0,0,0.2),(1,0,0,0.3)],filled=True,slices=[0])
#F.show_contour('../w51/H2CO_11_C_C_tausummed_42to61.fits',levels=[2.0,5.5],colors=[(1,0,0,0.4),(1,0,0,0.6)],filled=True,slices=[0])
#F.show_contour('../w51/H2CO_11_C_C_tausummed_63to67.fits',levels=[2.5,10.],colors=[(0,0,1,0.4),(0,0,1,0.6)],filled=True,slices=[0])
#F.show_contour('../w51_c_c/H2CO_11_C_C_tausummed_68to71.fits',levels=[0.6,5.0],colors=[(1,0,0,0.4),(1,0,0,0.6)],filled=True,slices=[0])
F.show_contour(dpath+'H2CO_22_Ku_D_tausummed_52to58.fits',levels=[0.5,2,15.0],colors=[(0.6,0.1,0.6,0.3),(0.6,0.1,0.6,0.5)],filled=True,slices=[0])
F.show_contour(dpath+'H2CO_22_Ku_D_tausummed_63to67.fits',levels=[0.9,  15.0],colors=[(0.0,0.5,1.0,0.3),(0.0,0.5,1.0,0.5)],filled=True,slices=[0])
F.show_contour(dpath+'H2CO_22_Ku_D_tausummed_67to70.fits',levels=[0.4,  15.0],colors=[(1.0,0.5,0.0,0.3),(1.0,0.5,0.0,0.5)],filled=True,slices=[0])
#F.show_contour('../w51_c_c/H2CO_11_C_C_absorbsummed_42to61.fits',levels=[-0.2,-0.3,-0.4],colors=[(1,0,0,0.1),(1,0,0,0.2),(1,0,0,0.3)],filled=True,slices=[0])
#F.show_contour('../w51_c_c/H2CO_11_C_C_absorbsummed_63to67.fits',levels=[-0.2,-0.3,-0.4],colors=[(0,1,1,0.1),(0,1,1,0.2),(0,1,1,0.3)],filled=True,slices=[0])
#F.show_contour('../w51_c_c/H2CO_11_C_C_absorbsummed_68to71.fits',levels=[-0.2,-0.3,-0.4],colors=[(1,0,1,0.1),(1,0,1,0.2),(1,0,1,0.3)],filled=True,slices=[0])
F.show_contour(dpath+'v2.0_ds2_l050_13pca_map20.fits',levels=[12],colors=['b'], linewidths=[2], filled=False,convention='calabretta')

#F.show_contour(dpath+'SgrB2_nh3_3-3_maximum.fits',levels=[0.04],colors=['b'])
#F.show_contour(dpath+'SGRB2_1.3CM_fix_gal.fits',levels=[0.003],colors=['r'],smooth=3)
#F.show_rectangles(sgrb2m.l.value-0.01,sgrb2m.b.value-0.005,width=5/60.,height=5/60.,color='g')

#F.show_regions(dpath+'scalebars_8.5kpc_gal.reg')
#F.show_regions(dpath+'nh3_observed_region.reg')

F.add_scalebar(length=((0.5 * u.pc)/(5.4*u.kpc)*u.radian).to(u.degree).value)
F.scalebar.set_label('0.5 pc')
F.scalebar.set_color('orange')
F.scalebar.set_linewidth(3)
F.scalebar.set_font_size(20)

beamd=(u.radian*((220*u.GHz).to(u.m,u.spectral()))/(12*u.m)).to(u.degree)
print beamd.to(u.arcsec)
beamd=beamd.value
F.add_beam(major=beamd/2,minor=beamd/2,hatch='///',facecolor='none',color=(1,1,0,1))
F.add_beam(major=0.5/3600.,minor=0.5/3600.,hatch='|||',facecolor='none',color=(0,1,1,1))

#F.show_regions('/Users/adam/work/w51/HCHII_candidates.reg')
F.show_regions('/Users/adam/work/w51/cycle2_ALMA_frame2.reg')

(xl,yl),(xu,yu) = F._ax1.bbox._bbox.corners()[[0,3]]

inset = aplpy.FITSFigure(hdu,convention='calabretta',figure=figure,subplot=[xl,yu-0.25,0.2,0.25])
inset.show_grayscale(stretch='arcsinh',vmin=-5e-4,vmax=0.031)
e1e2 = coordinates.ICRS(290.93268,14.508363,unit=('deg','deg'))
inset.recenter(e1e2.ra.value,e1e2.dec.value,width=15/60./60.,height=15/60./60.)
inset.tick_labels.hide()
inset.axis_labels.hide()
inset.show_contour(dpath+'W51Ku_BD_spw19.bigish_uniform_contsub19.clean.image.integ_52to65.fits', levels=[0.005,0.01], colors=['r','r'], smooth=3, slices=[0], linewidths=[2,2])
inset.add_beam(major=0.5/3600., minor=0.5/3600., color='orange', linewidth=3, label='1"')
inset.add_scalebar(length=((0.1 * u.pc)/(5.4*u.kpc)*u.radian).to(u.degree).value, color='orange', linewidth=3, label='0.1 pc')


F.save('/Users/adam/proposals/alma/cycle2/w51/W51_Ku_withH2COcontours.png')
F.save('/Users/adam/proposals/alma/cycle2/w51/W51_Ku_withH2COcontours.pdf')
