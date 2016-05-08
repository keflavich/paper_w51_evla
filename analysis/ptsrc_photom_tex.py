from __future__ import print_function
import numpy as np
from astropy import units as u
from astropy import table
import paths
from astropy import time
from rounded import rounded
from latex_info import latexdict, format_float
from astropy.table import Table
import copy

latexdict = copy.copy(latexdict)
tbl = Table.read(paths.tpath('EVLA_VLA_PointSourcePhotometry.ipac'), format='ascii.ipac')

## Create latex table

tbl.add_column(table.Column(data=tbl['peak_flux']-tbl['cutout_min_flux'],
                            name='peak_m_background',
                            unit=tbl['peak_flux'].unit))
nondetections = tbl['peak_m_background'] < tbl['local_rms_noise']*3

cols_order = ['SourceName', 'Epoch', 'ObservationDate', 'peak_flux', 'peak_m_background', 'local_rms_noise', 'Frequency']
cols = {'SourceName': 'Object',
        #'FrequencyName': 'Band',
        'ObservationDate': 'Obs. Date',
        'peak_flux': 'Peak $S_{\\nu}$',
        'peak_m_background': 'Peak - Background',
        'local_rms_noise': 'RMS',
       }

textbl = tbl.copy()[cols_order]
textbl.sort(['SourceName', 'Frequency'])
textbl[nondetections]['peak_flux'] = np.nan
textbl[nondetections]['peak_m_background'] = np.nan
textbl['peak_flux'] = ((list(map(lambda x,y: rounded(x,y,extra=0)[0],
                                 textbl['peak_flux'].to(u.mJy/u.beam).value,
                                 textbl['local_rms_noise'].to(u.mJy/u.beam).value))))
textbl['peak_m_background'] = ((list(map(lambda x,y: rounded(x,y,extra=0)[0],
                                         textbl['peak_m_background'].to(u.mJy/u.beam).value,
                                         textbl['local_rms_noise'].to(u.mJy/u.beam).value))))
textbl['local_rms_noise'] = ((list(map(lambda x,y: rounded(x,y,extra=0)[0],
                                       textbl['local_rms_noise'].to(u.mJy/u.beam).value,
                                       textbl['local_rms_noise'].to(u.mJy/u.beam).value))))
for name in ('peak_flux', 'peak_m_background', 'local_rms_noise'):
    textbl[name].unit = u.mJy/u.beam
textbl['SourceName'] = list(map(lambda x: x.replace("_","-"), textbl['SourceName']))

for old,new in cols.items():
    textbl.rename_column(old, new)

latexdict['header_start'] = '\label{tab:contsrcs}'
latexdict['caption'] = 'Continuum Point Sources (excerpt)'
latexdict['tablefoot'] = ('\par\nAn excerpt from the point source catalog.  '
                          'For the full catalog, see Table \\ref{tab:contsrcs_full}'
                          '')
print("latexdict for non-full version: ",latexdict)
textbl[::10].write(paths.tpath('pointsource_photometry.tex'),
                   format='ascii.latex',
                   latexdict=latexdict,
                   formats={'RMS': format_float,
                            'Peak $S_{\\nu}$': format_float,
                            'Peak - Background': format_float,
                            'Obs. Date': lambda x: time.Time(x).iso[:10],
                           })

latexdict['header_start'] = '\label{tab:contsrcs_full}'
latexdict['caption'] = 'Continuum Point Sources'
latexdict['tablefoot'] = ''
latexdict['tabletype'] = 'longtable'
latexdict['tablealign'] = 'c' * len(textbl.columns)
latexdict['col_align'] = ""
#latexdict['tabulartype'] = 'longtable'
from astropy.io.ascii.latex import Latex, LatexData, LatexHeader
from astropy.io import ascii

# allows overwriting the header type...
class LongTableData(LatexData):
    pass
class LongTableHeader(LatexHeader):
    pass
class LongTable(Latex):
    header_class = LongTableHeader
    data_class = LongTableData

LongTable.data_class.data_start = None
LongTable.data_class.data_end = "" #r'\end{longtable}'
LongTable.header_class.header_start = r"\\" #r'\begin{longtable}'

# unfortunately this table still requries some manipulation after the fact
# a \\ needs to be added after the \caption{}, since longtable doesn't seem to respect the caption command
ascii.write(textbl, output=paths.tpath('pointsource_photometry_full.tex'),
#ascii.write(textbl, output=paths.tpath('test.tex'),
            Writer=LongTable,
            #format='latex',
            latexdict=latexdict,
            formats={'RMS': format_float,
                     'Peak $S_{\\nu}$': format_float,
                     'Peak - Background': format_float,
                     'Obs. Date': lambda x: time.Time(x).iso[:10],
                    })

with open(paths.tpath('pointsource_photometry_full.tex'),'r') as fh:
    lines = fh.readlines()

with open(paths.tpath('pointsource_photometry_full.tex'),'w') as fh:
    lines[:3] = [r'\begin{longtable}{ccccccc}'+"\n",
                 r'\caption{Continuum Point Sources}\\'+"\n",
                 r''+"\n",]
    fh.writelines(lines)
