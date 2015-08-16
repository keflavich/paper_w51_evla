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

epoch3 = fits.open(paths.dpath("W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits"))
beam3 = Beam.from_fits_header(epoch3[0].header)
epoch3[0].data = epoch3[0].data.squeeze()
wcs3 = wcs.WCS(epoch3[0].header).sub([wcs.WCSSUB_CELESTIAL])
epoch3header = wcs3.to_header()
epoch3header['NAXIS'] = 2
epoch3header['NAXIS1'] = epoch3[0].data.shape[1]
epoch3header['NAXIS2'] = epoch3[0].data.shape[0]

kernel = (0.3*u.arcsec).to(u.deg)
pixscale = (wcs3.pixel_scale_matrix.diagonal()**2).sum()**0.5 * u.deg
print('pix kernel size: {0}'.format(kernel/pixscale))
smoothed3 = convolution.convolve_fft(epoch3[0].data,
                                     convolution.Gaussian2DKernel(kernel/pixscale))

diff3sm =  epoch3[0].data - smoothed3
diff3smhdu = fits.PrimaryHDU(data=diff3sm, header=epoch3header)
diff3smhdu.writeto(paths.dpath("Kuband_Epoch3sm-Epoch3.fits"), clobber=True, output_verify='fix')
