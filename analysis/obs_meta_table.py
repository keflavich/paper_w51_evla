import numpy as np
from astropy.table import Table
import paths
from latex_info import latexdict, exp_to_tex, format_float
import natsort

tbl = Table.read('../tables/obs_meta.tbl', format='ascii', delimiter="|")

tbl = tbl[natsort.index_natsorted(tbl['Date'])]

latexdict['header_start'] = '\label{tab:obs_meta}'
latexdict['caption'] = 'Observation Metadata'
latexdict['tablefoot'] = ('')
latexdict['col_align'] = 'l'*8
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
tbl.write(paths.tpath('obs_meta.tex'), format='ascii.latex',
          latexdict=latexdict, formats={},)

