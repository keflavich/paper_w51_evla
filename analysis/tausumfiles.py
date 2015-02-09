from paths import dpath,rpath,fpath
from astropy.io import fits

tau2 = fits.getdata(dpath('H2CO_22_Ku_D_taucube.fits'))
f2 = fits.open(dpath('H2CO_22_Ku_D_taucube.fits'))

# TODO below here
integ3 = tau2[0:62,:,:].sum(axis=0)
f2[0].data = integ3
f2[0].writeto('H2CO_22_Ku_D_tausummed_52to70.fits',clobber=True)
dpath = '/Volumes/128gbdisk/w51_c_c/'
f = fits.open(dpath+'H2CO_11_C_C_spw17_taucube.fits')
tau = f[0].data
integ3 = tau[18:62,:,:].sum(axis=0)
f[0].data = integ3
f[0].writeto('H2CO_11_C_C_tausummed_42to67.fits')
f = fits.open(dpath+'H2CO_11_speccube_contsub17_big_uniform.image.fits')
cspec = f[0].data[:70,:,:]
integ3 = cspec[18:62,:,:].sum(axis=0)
f[0].data = integ3
f[0].writeto('H2CO_11_C_C_absorbsummed_42to67.fits')
get_ipython().magic(u'pwd ')
imshow(tau2[30,:,:])
from pylab import *
imshow(tau2[30,:,:])
imshow(tau2[31,:,:])
imshow(tau2[35,:,:])
imshow(tau2[40,:,:])
imshow(np.isnan(tau2[40,:,:]))
imshow(np.isnan(tau2[20,:,:]))
imshow(np.isnan(tau2[21,:,:]))
np.isnan(tau2).sum(axis=1).sum(axis=1)
imshow(tau2[-9,:,:])
imshow(tau2[-10,:,:])
clf(); imshow(tau2[-10,:,:]); colorbar()
clf(); imshow(tau2[-9,:,:]); colorbar()
clf(); imshow(tau2[-11,:,:]); colorbar()
clf(); imshow(tau2[-8,:,:]); colorbar()
tau2.shape
clf(); imshow(tau2[50,:,:]); colorbar()
clf(); imshow(tau2[55,:,:]); colorbar()
clf(); imshow(tau2[56,:,:]); colorbar()
clf(); imshow(tau2[57,:,:]); colorbar()
clf(); imshow(tau2[58,:,:]); colorbar()
clf(); imshow(tau2[54,:,:]); colorbar()
clf(); imshow(tau2[53,:,:]); colorbar()
tau2[np.isnan(tau2)] = 9
clf(); imshow(tau2[53,:,:]); colorbar()
tau2 = fits.getdata('H2CO_22_Ku_D_taucube.fits')
tau2[np.isnan(tau2)] = 5
clf(); imshow(tau2[53,:,:]); colorbar()
tau2 = fits.getdata('H2CO_22_Ku_D_taucube.fits')
tau2[np.isnan(tau2)] = 1
clf(); imshow(tau2[53,:,:]); colorbar()
tau2 = fits.getdata('H2CO_22_Ku_D_taucube.fits')
tau2[np.isnan(tau2)] = 1
integ3 = tau2[0:62,:,:].sum(axis=0)
f2[0].data = integ3
f2[0].writeto('H2CO_22_Ku_D_tausummed_52to70.fits',clobber=True)
get_ipython().magic(u'ls -rt')
yy,xx = np.indices(tau2.shape[1:])
yty
yy
yy,xx = np.indices(tau2.shape[1:])
rr = ((xx-255.)**2+(yy-255.)**2)**0.5
get_ipython().magic(u'paste')
yy,xx = np.indices(tau2.shape[1:])
rr = ((xx-255.)**2+(yy-255.)**2)**0.5
mask = rr < 178
nanmask = np.zeros_like(mask)
nanmask[True-mask] = np.nan

integ1 = tau2[0:23,:,:].sum(axis=0)
f2[0].data = integ1 + nanmask
f2[0].writeto('H2CO_22_Ku_D_tausummed_52to58.fits',clobber=True)

integ2 = tau2[38:51,:,:].sum(axis=0)
f2[0].data = integ2 + nanmask
f2[0].writeto('H2CO_22_Ku_D_tausummed_63to67.fits',clobber=True)

integ3 = tau2[51:62,:,:].sum(axis=0)
f2[0].data = integ3 + nanmask
f2[0].writeto('H2CO_22_Ku_D_tausummed_67to70.fits',clobber=True)


integ3 = tau2[0:62,:,:].sum(axis=0)
f2[0].data = integ3 + nanmask
f2[0].writeto('H2CO_22_Ku_D_tausummed_52to70.fits',clobber=True)
nanmask
nanmask = np.zeros_like(mask,dtype='float')
nanmask[True-mask] = np.nan
nanmask
get_ipython().magic(u'paste')
nanmask = np.zeros_like(mask,dtype='float')
nanmask[True-mask] = np.nan

integ1 = tau2[0:23,:,:].sum(axis=0)
f2[0].data = integ1 + nanmask
f2[0].writeto('H2CO_22_Ku_D_tausummed_52to58.fits',clobber=True)

integ2 = tau2[38:51,:,:].sum(axis=0)
f2[0].data = integ2 + nanmask
f2[0].writeto('H2CO_22_Ku_D_tausummed_63to67.fits',clobber=True)

integ3 = tau2[51:62,:,:].sum(axis=0)
f2[0].data = integ3 + nanmask
f2[0].writeto('H2CO_22_Ku_D_tausummed_67to70.fits',clobber=True)


integ3 = tau2[0:62,:,:].sum(axis=0)
f2[0].data = integ3 + nanmask
f2[0].writeto('H2CO_22_Ku_D_tausummed_52to70.fits',clobber=True)
tau.sshape
tau.shape
########################################################
# Started Logging At: 2013-12-07 01:32:43
########################################################

########################################################
# # Started Logging At: 2013-12-07 01:32:44
########################################################
get_ipython().magic(u'run -i ~/repos/w51evlareductionscripts/make_taucube_vla_ku_bd_22.py')
imshow(cont2)
from pylab import *
imshow(cont2)
clf(); imshow(cont2); colorbar()
clf(); imshow(cont2,vmax=0.01); colorbar()
clf(); imshow(cont2,vmax=0.005); colorbar()
clf(); imshow(cont2,vmax=0.005,cmap=matplotlib.cm.spectral); colorbar()
clf(); imshow(cont2,vmax=0.001,cmap=matplotlib.cm.spectral); colorbar()
clf(); imshow(cont2==0,cmap=matplotlib.cm.spectral); colorbar()
clf(); imshow(cont2==cont2.min(),cmap=matplotlib.cm.spectral); colorbar()
cont2.min()
d=fits.getdata(dpath+'W51Ku_BD_spw20.small_uniform_continuum.clean.image.fits').squeeze()
d[100:200,300:500].std()
get_ipython().magic(u'run -i ~/repos/w51evlareductionscripts/make_taucube_vla_ku_bd_22.py')
import agpy
agpy.smooth
get_ipython().magic(u'run -i ~/repos/w51evlareductionscripts/make_taucube_vla_ku_bd_22.py')
get_ipython().magic(u'pinfo smooth')
get_ipython().magic(u'pinfo agpy.smooth')
get_ipython().magic(u'run -i ~/repos/w51evlareductionscripts/make_taucube_vla_ku_bd_22.py')
tau2 = -np.log((spec2s/cont2))
tau2[np.isnan(tau2)] = 5
f2[0].header['CRPIX3'] -= 190
f2[0].data = tau2
f2.writeto('H2CO_22_Ku_BD_small_taucube_sm.fits',clobber=True)
get_ipython().magic(u'run -i ~/repos/w51evlareductionscripts/make_taucube_vla_ku_bd_22.py')
get_ipython().magic(u'paste')
tau2 = -np.log((spec2s/cont2s))
tau2[np.isnan(tau2)] = 5

f2[0].header['CRPIX3'] -= 190
f2[0].data = tau2
f2.writeto('H2CO_22_Ku_BD_small_taucube_sm.fits',clobber=True)
cont2s+=0.05
get_ipython().magic(u'paste')
#tau2 = -np.log((spec2/cont2))
tau2 = -np.log((spec2s/cont2s))
tau2[np.isnan(tau2)] = 5

f2[0].header['CRPIX3'] -= 190
f2[0].data = tau2
f2.writeto('H2CO_22_Ku_BD_small_taucube_sm.fits',clobber=True)
spec2s += contlev+0.05
#tau2 = -np.log((spec2/cont2))
tau2 = -np.log((spec2s/cont2s))
tau2[np.isnan(tau2)] = 5
f2[0].header['CRPIX3'] -= 190
f2[0].data = tau2
f2.writeto('H2CO_22_Ku_BD_small_taucube_sm.fits',clobber=True)
get_ipython().magic(u'run -i ~/repos/w51evlareductionscripts/make_taucube_vla_ku_bd_22.py')
cont2s = agpy.smooth(cont2+contlev,kernelwidth=2)
tau2 = -np.log((spec2s/cont2s))
tau2[np.isnan(tau2)] = 5
f2[0].header['CRPIX3'] -= 190
f2[0].data = tau2
f2.writeto('H2CO_22_Ku_BD_small_taucube_sm.fits',clobber=True)
get_ipython().magic(u'run -i ~/repos/w51evlareductionscripts/make_taucube_vla_ku_bd_22.py')
get_ipython().magic(u'paste')
yy,xx = np.indices(tau2.shape[1:])
rr = ((xx-255.)**2+(yy-255.)**2)**0.5
mask = rr < 200
nanmask = np.zeros_like(mask,dtype='float')
nanmask[True-mask] = np.nan

integ1 = tau2[11:48,:,:].sum(axis=0)
f2[0].data = integ1 + nanmask
f2[0].writeto('H2CO_22_Ku_BD_small_tausummed_48to61.fits',clobber=True)

integ2 = tau2[42:70,:,:].sum(axis=0)
f2[0].data = integ2 + nanmask
f2[0].writeto('H2CO_22_Ku_BD_small_tausummed_59to68.fits',clobber=True)

integ3 = tau2[70:81,:,:].sum(axis=0)
f2[0].data = integ3 + nanmask
f2[0].writeto('H2CO_22_Ku_BD_small_tausummed_68to72.fits',clobber=True)


integ3 = tau2[11:81,:,:].sum(axis=0)
f2[0].data = integ3 + nanmask
f2[0].writeto('H2CO_22_Ku_BD_small_tausummed_48to72.fits',clobber=True)
81-11
tau2ds = (tau2[11:83::3,:,:]+tau2[12:83::3,:,:]+tau2[13:83::3,:,:])/3.
tau2ds = (tau2[11:83:3,:,:]+tau2[12:83:3,:,:]+tau2[13:83:3,:,:])/3.
h = f2[0].header
(11-h['CRPIX3']+1)*h['CDELT3']+h['CRVAL3']
(11-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']
get_ipython().magic(u'paste')
tau2ds = (tau2[11:83:3,:,:]+tau2[12:83:3,:,:]+tau2[13:83:3,:,:])/3.
h = f2[0].header
f2[0].header['CRVAL3'] = (10-h['CRPIX3']+1)*h['CDELT3']+h['CRVAL3']
f2[0].header['CRPIX3'] = 1
f2[0].header['CDELT3'] = h['CDELT3']*3
f2[0].data = tau2ds
f2.writeto('H2CO_22_Ku_BD_small_taucube_sm_ds.fits',clobber=True)

yy,xx = np.indices(tau2.shape[1:])
rr = ((xx-255.)**2+(yy-255.)**2)**0.5
mask = rr < 200
nanmask = np.zeros_like(mask,dtype='float')
nanmask[True-mask] = np.nan

integ1 = tau2[10:47,:,:].sum(axis=0)
f2[0].data = integ1 + nanmask
f2[0].writeto('H2CO_22_Ku_BD_small_tausummed_48to61.fits',clobber=True)

integ2 = tau2[41:69,:,:].sum(axis=0)
f2[0].data = integ2 + nanmask
f2[0].writeto('H2CO_22_Ku_BD_small_tausummed_59to68.fits',clobber=True)

integ3 = tau2[69:80,:,:].sum(axis=0)
f2[0].data = integ3 + nanmask
f2[0].writeto('H2CO_22_Ku_BD_small_tausummed_68to72.fits',clobber=True)


integ3 = tau2[10:80,:,:].sum(axis=0)
f2[0].data = integ3 + nanmask
f2[0].writeto('H2CO_22_Ku_BD_small_tausummed_48to72.fits',clobber=True)
tau2ds.shape
