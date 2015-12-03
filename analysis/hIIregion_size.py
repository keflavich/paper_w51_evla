import numpy as np
from astropy.io import fits
from astropy import units as u
import radio_beam
import paths
from common_constants import distance

header = fits.getheader(paths.dpath("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits"))
beam = radio_beam.Beam.from_fits_header(header)

mJytoK = (1*u.mJy).to(u.K, u.brightness_temperature(beam, 4.9*u.GHz))

print("A 1 mJy source will have T_B = {0} at 4.9 GHz".format(mJytoK))
sizelim = (((mJytoK/(1e4*u.K) * (beam.sr))/(2.*np.pi))**0.5).to(u.arcsec)
print("The upper size limit for a 10^4 K source is {0},"
      "or {1} = {2}".format(sizelim, (sizelim*distance).to(u.mpc,
                                                           u.dimensionless_angles()),
                            (sizelim*distance).to(u.au,
                                                  u.dimensionless_angles())
                           ))
