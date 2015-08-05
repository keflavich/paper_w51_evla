from astropy.io import fits
from astropy import coordinates
from astropy import units as u
from astropy import wcs
import FITS_tools
import reproject
import paths
import image_registration
from radio_beam import Beam
from astropy import convolution
import gaussfitter

epoch1 = fits.open(paths.dpath("W51-CBAND-feathered.fits"))
epoch1[0].data = epoch1[0].data.squeeze()
beam1 = Beam.from_fits_header(epoch1[0].header)
epoch1[0].header = epoch1header = wcs.WCS(epoch1[0].header).sub([wcs.WCSSUB_CELESTIAL]).to_header()
epoch1header['NAXIS1'] = epoch1[0].data.shape[1]
epoch1header['NAXIS2'] = epoch1[0].data.shape[0]
wcs1 = wcs.WCS(epoch1[0].header).sub([wcs.WCSSUB_CELESTIAL])
epoch3 = fits.open(paths.dpath("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits"))
beam3 = Beam.from_fits_header(epoch3[0].header)
epoch3[0].data = epoch3[0].data.squeeze()
wcs3 = wcs.WCS(epoch3[0].header).sub([wcs.WCSSUB_CELESTIAL])
epoch3header = wcs3.to_header()
epoch3header['NAXIS'] = 2
epoch3header['NAXIS1'] = epoch3[0].data.shape[1]
epoch3header['NAXIS2'] = epoch3[0].data.shape[0]

#epoch1reproj = FITS_tools.hcongrid.hcongrid(epoch1[0].data, epoch1header, epoch3header)
epoch1reproj, footprint = reproject.reproject(epoch1[0], epoch3header)

scalefactor = (beam3.sr/beam1.sr).value

xshift,yshift, ex, ey = image_registration.chi2_shift(epoch3[0].data, epoch1reproj*scalefactor, err=0.001)

center = coordinates.SkyCoord(290.92443*u.deg, 14.515755*u.deg, frame='fk5')
xx,yy = wcs1.wcs_world2pix([[center.fk4.ra.deg, center.fk4.dec.deg]], 0)[0]
subim1 = epoch1[0].data[yy-10:yy+10, xx-10:xx+10]
wcs1sub = wcs1[yy-10:yy+10, xx-10:xx+10]
xx,yy = wcs3.wcs_world2pix([[center.fk5.ra.deg, center.fk5.dec.deg]], 0)[0]
wcs3sub = wcs3[yy-10:yy+10, xx-10:xx+10]
subim3 = epoch3[0].data[yy-10:yy+10, xx-10:xx+10]
gf1 = gaussfitter.gaussfit(subim1)
gf3 = gaussfitter.gaussfit(subim3)

dx1,dy1 = gf1[2:4]
dx3,dy3 = gf3[2:4]
dx_gf, dy_gf = dx3-dx1, dy3-dy1

epoch1matched_fft = image_registration.fft_tools.shift2d(epoch1reproj*scalefactor, -xshift, -yshift)
epoch1matched = image_registration.fft_tools.shift2d(epoch1reproj*scalefactor, dx_gf, dy_gf)

diff =  epoch3[0].data - epoch1matched
diffhdu = fits.PrimaryHDU(data=diff, header=epoch3header)
diffhdu.writeto(paths.dpath("Cband_Epoch3-Epoch1.fits"), clobber=True, output_verify='fix')

smooth_beam = beam1.deconvolve(beam3)
kernel = smooth_beam.major
pixscale = (wcs3.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg
smoothed3 = convolution.convolve_fft(epoch3[0].data,
                                     convolution.Gaussian2DKernel(kernel/pixscale))

diffsm =  smoothed3 - epoch1matched
diffsmhdu = fits.PrimaryHDU(data=diffsm, header=epoch3header)
diffsmhdu.writeto(paths.dpath("Cband_Epoch3sm-Epoch1.fits"), clobber=True, output_verify='fix')

diff3sm =  epoch3[0].data - smoothed3
diff3smhdu = fits.PrimaryHDU(data=diff3sm, header=epoch3header)
diff3smhdu.writeto(paths.dpath("Cband_Epoch3sm-Epoch3.fits"), clobber=True, output_verify='fix')

import pylab as pl
import mpl_plot_templates
pl.figure(1).clf()
mpl_plot_templates.adaptive_param_plot(smoothed3.ravel(), epoch1matched.ravel())
pl.plot([0,0.012], [0,0.012], 'k--', alpha=0.3)
pl.figure(2).clf()
mpl_plot_templates.adaptive_param_plot(smoothed3.ravel(), epoch1matched_fft.ravel())
pl.plot([0,0.012], [0,0.012], 'k--', alpha=0.3)

