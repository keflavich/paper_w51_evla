import numpy as np
import paths
import pyregion
from astropy import coordinates
from astropy import units as u
from astropy import table
from astropy.table import Table,Column
import latex_info
reload(latex_info)
from latex_info import latexdict, exp_to_tex, format_float
import radio_beam
from astropy.io import fits

from photom_files import files

beams = {name: radio_beam.Beam.from_fits_header(fits.getheader(fn))
         for name,fn in files.items()}
dates = {name: fits.getheader(fn)['DATE-OBS']
         for name,fn in files.items()}

obstbl = Table([Column(data=[ep.split()[-1] for ep in beams], name='Epoch'),
                Column(data=[float(ep.split()[0]) for ep in beams]*u.GHz, name='Frequency'),
                Column(data=[bm.major.to(u.arcsec).value for bm in beams.values()]*u.arcsec, name='BMAJ'),
                Column(data=[bm.minor.to(u.arcsec).value for bm in beams.values()]*u.arcsec, name='BMIN'),
                Column(data=[bm.pa.to(u.deg).value for bm in beams.values()]*u.deg, name='BPA'),
               ])

obstbl.sort(['Frequency', 'Epoch'])
latexdict['header_start'] = '\label{tab:observations}'
latexdict['caption'] = 'Observations'
#latexdict['col_align'] = 'lllrr'
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
obstbl.write(paths.tpath('observations.tex'), format='ascii.latex',
             latexdict=latexdict,
             formats={'BMAJ': lambda x: '{0:0.2f}'.format(x),
                      'BMIN': lambda x: '{0:0.2f}'.format(x),
                      'BPA':  lambda x: '{0:0.2f}'.format(x),
                     },
            )
 
