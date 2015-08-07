import FITS_tools
from astropy.io import fits
import paths

neiir = FITS_tools.regrid_cube_hdu(fits.open(paths.dpath('w51.neii.fits'))[0], outheader=fits.Header.fromtextfile(paths.dpath('w51.neii.square.hdr')))
neiir.writeto(paths.dpath('w51.neii.square.fits'))

sivr = FITS_tools.regrid_cube_hdu(fits.open(paths.dpath('w51.siv.fits'))[0], outheader=fits.Header.fromtextfile(paths.dpath('w51.siv.square.hdr')))
sivr.writeto(paths.dpath('w51.siv.square.fits'))
