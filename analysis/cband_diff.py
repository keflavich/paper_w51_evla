from astropy.io import fits
from astropy import wcs
import FITS_tools
import reproject
import paths
import image_registration
from radio_beam import Beam
from astropy import convolution

epoch1 = fits.open(paths.dpath("W51-CBAND-feathered.fits"))
epoch1[0].data = epoch1[0].data.squeeze()
beam1 = Beam.from_fits_header(epoch1[0].header)
epoch1[0].header = epoch1header = wcs.WCS(epoch1[0].header).sub([wcs.WCSSUB_CELESTIAL]).to_header()
epoch1header['NAXIS1'] = epoch1[0].data.shape[1]
epoch1header['NAXIS2'] = epoch1[0].data.shape[0]
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

epoch1matched = image_registration.fft_tools.shift2d(epoch1reproj*scalefactor, -xshift, -yshift)

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
