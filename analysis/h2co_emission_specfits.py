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
                  "spectra/emission/W51Ku_BD_h2co_v30to90_natural_contsub.image*.fits")))
      if 'mol' not in x and '?' not in x
     ]
spectra = sp


tbl = Table(dtype=[('Object Name', 'str20'),
                   ('Amplitude', float),
                   ('$\sigma(Amplitude)', float),
                   ('$V_{LSR}$', float),
                   ('$\sigma(V_{LSR})$', float),
                   ('$dV$', float),
                   ('$\sigma(dV)', float),
                   ('$\Omega_{ap}$', float),]
           )

# My manual inspection: which are detected?
# weakdetections are those that are not clearly believable
detections = ['e8mol', 'e2-e8 bridge', 'e10mol', 'NorthCore']
weakdetections = []

# conversion....
[s.xarr.convert_to_unit('km/s') for s in sp]

# setup
for s in sp:
    s.specname = s.fileprefix.split("_")[-1]
    log.info(s.specname+" stats")
    noiseregion = (s.xarr < 20*u.km/u.s).value | (s.xarr > 100*u.km/u.s).value
    s.error[:] = s.data[noiseregion].std()

# fitting
for s in sp:
    s.plotter(xmin=-10,xmax=120)
    s.specfit(fittype='gaussian',
              guesses=[0.03,55,10],
              limited=[(True,False),(False,False),(True,False)])
    log.info(s.specname+" fitting: {0}".format(s.specfit.parinfo))
    s.plotter.ymin -= 0.005
    s.specfit.plotresiduals(axis=s.plotter.axis,clear=False,yoffset=-0.005,label=False)
    s.plotter.savefig(paths.fpath('spectra/emission/'+s.specname+"_h2co22emisson_fit.png"),
                                  bbox_inches='tight')

    tbl.add_row([sp.specname,
                 sp.specfit.parinfo.AMPLITUDE0.value, sp.specfit.parinfo.AMPLITUDE0.error,
                 sp.specfit.parinfo.SHIFT0.value, sp.specfit.parinfo.SHIFT0.error,
                 sp.specfit.parinfo.WIDTH0.value, sp.specfit.parinfo.WIDTH0.error,
                 sp.header['APAREA']])

 
# sort such that e10 comes after e9
import natsort
tbl = tbl[natsort.index_natsorted(tbl['Object Name'])]

detection_note = ['-' if name in detections else
                  'weak' if name in weakdetections else
                  'none'
                  for name in tbl['ObjectName']]
tbl.add_column(table.Column(data=detection_note, name='Detection Status'))

ok = np.array([row['ObjectName'] in detections+weakdetections
               for row in tbl])

tbl[ok].write(paths.tpath('H2CO22_emission_spectral_fits.ecsv'), format='ascii.ecsv')

tbl[ok].write(paths.tpath('H2CO22_emission_spectral_fits.tex'), format='ascii.latex')
