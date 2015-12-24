import numpy as np
import paths
import pyregion
from astropy import coordinates
from astropy import units as u
from astropy import table
from astropy.table import Table,Column
import latex_info
from latex_info import latexdict, exp_to_tex, format_float
import radio_beam
from astropy.io import fits
from astropy import wcs

from photom_files import files

# add images that merge subbands
files["5 GHz Epoch 3"] = paths.dpath("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits")
files["13 GHz Epoch 2"] = paths.dpath('W51Ku_BDarray_continuum_2048_both_uniform.hires.clean.image.fits')

beams = {name: radio_beam.Beam.from_fits_header(fits.getheader(fn))
         for name,fn in files.items()}
dates = {name: fits.getheader(fn)['DATE-OBS']
         for name,fn in files.items()}

noisereg = pyregion.open(paths.rpath("noise_estimate_region.reg"))

noise_est = {}
peak = {}

# noise estimate
for ep in files:
    im = fits.open(files[ep])
    header = im[0].header
    w = wcs.WCS(header).sub([wcs.WCSSUB_CELESTIAL])
    mask = noisereg.get_mask(header=w.to_header(), shape=im[0].data.squeeze().shape)
    noise_est[ep] = im[0].data.squeeze()[mask].std() * 1000
    peak[ep] = im[0].data.max()*1000
    

obstbl = Table([Column(data=[ep.split()[-1] for ep in beams], name='Epoch'),
                Column(data=[float(ep.split()[0]) for ep in beams]*u.GHz, name='Frequency'),
                Column(data=[bm.major.to(u.arcsec).value for bm in beams.values()]*u.arcsec, name='BMAJ'),
                Column(data=[bm.minor.to(u.arcsec).value for bm in beams.values()]*u.arcsec, name='BMIN'),
                Column(data=[bm.pa.to(u.deg).value for bm in beams.values()]*u.deg, name='BPA'),
                Column(data=[noise_est[ep] for ep in beams]*u.mJy, name='Noise Estimate'),
                Column(data=[peak[ep]/noise_est[ep] for ep in beams], name='Dynamic Range'),
               ])
obstbl.add_column(Column(data=[(1*u.Jy).to(u.K,
                                           u.brightness_temperature(bm.sr, freq*u.GHz,)
                                          ).value
                               for bm,freq in zip(beams.values(),
                                                  obstbl['Frequency'])
                              ],
                         name='Jy-Kelvin'
                        )
                 )

obstbl.sort(['Frequency', 'Epoch'])
latexdict['header_start'] = '\label{tab:observations}'
latexdict['caption'] = 'Observations'
latexdict['tablefoot'] = ('\par\nJy-Kelvin gives the conversion factor from Jy '
                          'to Kelvin given the synthesized beam size and '
                          'observation frequency')
#latexdict['col_align'] = 'lllrr'
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
obstbl.write(paths.tpath('observations.tex'), format='ascii.latex',
             latexdict=latexdict,
             formats={'BMAJ': lambda x: '{0:0.2f}'.format(x),
                      'BMIN': lambda x: '{0:0.2f}'.format(x),
                      'BPA':  lambda x: '{0:0.2f}'.format(x),
                      'Noise Estimate': lambda x: '{0:0.2f}'.format(x),
                      'Dynamic Range': lambda x: '{0:d}'.format(int(round(x))),
                      'Jy-Kelvin':  format_float,
                     },
            )
 
