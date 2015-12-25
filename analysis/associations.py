import numpy as np
import paths
from astropy import coordinates
from astropy import units as u
from astropy.table import Table,Column
from latex_info import latexdict, exp_to_tex, format_float

with open(paths.tpath("associations.tbl")) as f:
    lines = f.readlines()
    split1 = lines[0].find('Xray_Association')
    split2 = lines[0].find('NIR_Association')
    xrayassociations = {x[:split1].strip():x[split1:split2].strip() for x in lines[1:]}
    nirassociations = {x[:split1].strip():x[split2:].strip() for x in lines[1:]}
    names = sorted(xrayassociations.keys())
    keep = [name for name in names if not (xrayassociations[name]=='-' and nirassociations=='-')]
    xrayassociations_column = Column(data=[xrayassociations[x] for x in keep],
                                     name="X-ray")
    nirassociations_column = Column(data=[nirassociations[x] for x in keep],
                                    name="NIR")
    names_column = Column(data=keep, name='Source Name')

tbl = Table([names_column, xrayassociations_column, nirassociations_column])

import natsort
tbl = tbl[natsort.index_natsorted(tbl['Source Name'])]

latexdict['header_start'] = '\label{tab:associations}'
latexdict['caption'] = 'Source Associations'
latexdict['tablefoot'] = ('')
latexdict['col_align'] = 'lll'
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
tbl.write(paths.tpath('associations.tex'), format='ascii.latex',
          latexdict=latexdict, formats={},)

