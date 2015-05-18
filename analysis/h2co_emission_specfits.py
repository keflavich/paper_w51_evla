""" was vdisp.py """
import numpy as np
import pyspeckit 
import glob
from astropy import table
from astropy import units as u
from astropy import log
from astropy.utils.console import ProgressBar
import paths
from astropy.table import Table, Column

sp = [pyspeckit.Spectrum(x) for x in
      ProgressBar(
          glob.glob(
              paths.dpath(
                  "spectra/emission/W51Ku_BD_h2co_v30to90_briggs0_contsub.image*.fits")))
      if '?' not in x
     ]
spectra = sp


tbl = Table(dtype=[(str, 20), float,   float,   float,   float,   float,
                   float,   float],
                   names=['Object Name',
                          'Amplitude',
                          '$\sigma(Amplitude)',
                          '$V_{LSR}$',
                          '$\sigma(V_{LSR})$',
                          '$dV$',
                          '$\sigma(dV)',
                          '$\Omega_{ap}$',],
           )

# My manual inspection: which are detected?
# weakdetections are those that are not clearly believable
detections = ['e8mol', 'e2-e8 bridge', 'e10mol', 'NorthCore']
weakdetections = ['e8mol_ext', 'e10mol_ext']

# conversion....
[s.xarr.convert_to_unit('km/s') for s in sp]

# setup
for s in sp:
    s.specname = s.fileprefix.split("_")[-1]
    log.info(s.specname+" stats")
    noiseregion = (s.xarr < 40*u.km/u.s).value | (s.xarr > 80*u.km/u.s).value
    assert np.any(noiseregion)
    s.error[:] = s.data[noiseregion].std()

# fitting
for thisspec in sp:
    thisspec.plotter(xmin=30,xmax=90)
    thisspec.specfit(fittype='gaussian',
              guesses=[0.03,55,3],
              limited=[(True,False),(False,False),(True,False)])
    log.info(thisspec.specname+" fitting: {0}".format(thisspec.specfit.parinfo))
    thisspec.plotter.ymin -= 0.005
    thisspec.specfit.plotresiduals(axis=thisspec.plotter.axis,clear=False,yoffset=-0.005,label=False)
    thisspec.plotter.savefig(paths.fpath('spectra/emission/'+thisspec.specname+"_h2co22emisson_fit.png"),
                                  bbox_inches='tight')

    tbl.add_row([thisspec.specname,
                 thisspec.specfit.parinfo.AMPLITUDE0.value, thisspec.specfit.parinfo.AMPLITUDE0.error,
                 thisspec.specfit.parinfo.SHIFT0.value, thisspec.specfit.parinfo.SHIFT0.error,
                 thisspec.specfit.parinfo.WIDTH0.value, thisspec.specfit.parinfo.WIDTH0.error,
                 thisspec.header['APAREA']])

 
# sort such that e10 comes after e9
import natsort
tbl = tbl[natsort.index_natsorted(tbl['Object Name'])]

detection_note = ['-' if name in detections else
                  'weak' if name in weakdetections else
                  'none'
                  for name in tbl['Object Name']]
tbl.add_column(table.Column(data=detection_note, name='Detection Status'))

ok = np.array([row['Object Name'] in detections+weakdetections
               for row in tbl])

tbl[ok].write(paths.tpath('H2CO22_emission_spectral_fits.ecsv'), format='ascii.ecsv')

tbl[ok].write(paths.tpath('H2CO22_emission_spectral_fits.tex'), format='ascii.latex')
