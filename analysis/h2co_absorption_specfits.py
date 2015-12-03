""" was vdisp.py """
import numpy as np
import pyspeckit 
import glob
from astropy.io import ascii
from astropy import table
from astropy import units as u
from astropy import log
from astropy.utils.console import ProgressBar
import paths
from astropy.table import Table, Column
from rounded import rounded
from latex_info import latexdict, format_float
import pylab as pl

sp = [pyspeckit.Spectrum(x) for x in
      ProgressBar(
          glob.glob(
              paths.dpath(
                  "spectra/hiiregionh2co/W51Ku_BD_h2co_v30to90_briggs0_contsub.image*.fits")))
      if '?' not in x and 'mol' not in x
     ]
spectra = sp


tbl = Table(dtype=[(str, 20), float,   float,   float,   float,   float,
                   float,   float],
                   names=['Object Name',
                          'Amplitude',
                          '$E$(Amplitude)',
                          '$V_{LSR}$',
                          '$E(V_{LSR})$',
                          '$\sigma_V$',
                          '$E(\sigma_V)$',
                          '$\Omega_{ap}$',],
           )

# My manual inspection: which are detected?
# ambiguous are those that are not clearly believable
detections = ['e2', 'e3', 'e5', 'e6', 'e9', 'e10']
ambiguousdetections = ['e1']#'e8mol_ext', 'e10mol_ext']

# conversion....
[s.xarr.convert_to_unit('km/s') for s in sp]

# setup
for s in sp:
    assert s.specname
    log.info(s.specname+" stats")
    noiseregion = (s.xarr < 40*u.km/u.s).value | (s.xarr > 80*u.km/u.s).value
    assert np.any(noiseregion)
    s.error[:] = s.data[noiseregion].std()

# fitting
for thisspec in sp:
    fig = pl.figure(1)
    fig.clf()
    thisspec.plotter(xmin=30,xmax=90, figure=fig)
    thisspec.specfit(fittype='gaussian',
              guesses=[-0.03,thisspec.xarr[thisspec.data.argmin()].value,3],
              limited=[(False,True),(False,False),(True,False)])
    log.info(thisspec.specname+" fitting: {0}".format(thisspec.specfit.parinfo))
    thisspec.plotter.ymin -= 0.005
    thisspec.plotter.ymax += 0.005
    thisspec.specfit.annotate(loc='lower left')
    thisspec.specfit.plotresiduals(axis=thisspec.plotter.axis,clear=False,yoffset=+0.005,label=False)
    thisspec.plotter.savefig(paths.fpath('spectra/hiiregionh2co/'+thisspec.specname+"_h2co22absorption_fit.png"),
                                  bbox_inches='tight')

    tbl.add_row([thisspec.specname,]+
                 list((rounded(thisspec.specfit.parinfo.AMPLITUDE0.value, thisspec.specfit.parinfo.AMPLITUDE0.error)*u.Jy).to(u.mJy))+
                 list(rounded(thisspec.specfit.parinfo.SHIFT0.value, thisspec.specfit.parinfo.SHIFT0.error)*u.km/u.s)+
                 list(rounded(thisspec.specfit.parinfo.WIDTH0.value, thisspec.specfit.parinfo.WIDTH0.error)*u.km/u.s)+
                 [np.round(thisspec.header['APAREA']*(np.pi/180.)**2, int(np.ceil(-np.log10(thisspec.header['APAREA']*(np.pi/180.)**2)))+1)*u.sr])

 
# sort such that e10 comes after e9
import natsort
tbl = tbl[natsort.index_natsorted(tbl['Object Name'])]

detection_note = ['-' if name in detections else
                  'ambig' if name in ambiguousdetections else
                  'none'
                  for name in tbl['Object Name']]
tbl.add_column(table.Column(data=detection_note, name='Detection Status'))

ok = np.array([row['Object Name'] in detections+ambiguousdetections
               for row in tbl])

tbl[ok].write(paths.tpath('H2CO22_hiiregion_spectral_fits.ecsv'), format='ascii.ecsv')

for row in tbl:
    if "_" in row['Object Name']:
        row['Object Name'] = row['Object Name'].replace("_","\_")

latexdict['header_start'] = '\label{tab:absorption22}'
latexdict['caption'] = '\\formaldehyde \\twotwo absorption line parameters'
tbl[ok].write(paths.tpath('H2CO22_hiiregion_spectral_fits.tex'), format='ascii.latex', latexdict=latexdict,
              formats={'$\Omega_{ap}$': format_float}
             )
