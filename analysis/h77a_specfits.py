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
from rounded import rounded
from latex_info import latexdict

# Read in all spectra extracted to the h77a directory matching the "best" h77a file
# May 18, 2015: I think H77a big2 is best, but I haven't inspected.  EDIT: big2 is in velocity.
# Exclude "mol" files since they are molecular and should be extracted from the molecular cubes
# Exclude "?" files since they are misnamed and use the obsolete scheme from
# before I decided upon e9/e10
sp = [pyspeckit.Spectrum(x) for x in
      ProgressBar(
          glob.glob(
              paths.dpath(
                  "spectra_h77/H77a_BDarray_speccube_uniform_contsub_cvel_big2*e[0-9]*.fits")))
      if 'mol' not in x and '?' not in x
     ]
spectra = sp

# My manual inspection: which are detected?
# weakdetections are those that are not clearly believable
detections = ['e1', 'e2', 'e3', 'e4', 'e6',]
weakdetections = ['e5', 'e9', 'e10']

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
              guesses=[1,55,10],
              limited=[(True,False),(False,False),(True,False)])
    log.info(s.specname+" fitting: {0}".format(s.specfit.parinfo))
    s.plotter.ymin -= 0.0003
    s.specfit.plotresiduals(axis=s.plotter.axis,clear=False,yoffset=-0.0003,label=False)
    s.plotter.savefig(paths.fpath('spectra/h77/'+s.specname+"_h77a_fit.png"),
                                  bbox_inches='tight')

tbl = table.Table()
names = table.Column(data=[sp.specname for sp in spectra], name='ObjectName')
tbl.add_column(names)
for ii,(parname,unit) in enumerate([('amplitude',u.mJy/u.beam),
                             ('velocity',u.km/u.s),
                             ('width',u.km/u.s)]):
    dataerror = [rounded(sp.specfit.parinfo[ii].value, sp.specfit.parinfo[ii].error)
                 for sp in spectra]
    data,error = zip(*dataerror)
    #data = [sp.specfit.parinfo[ii].value
    #        for sp in spectra]
    #error = [sp.specfit.parinfo[ii].error
    #        for sp in spectra]
    column = table.Column(data=data,
                          name='H77a_'+parname,
                          unit=unit)
    tbl.add_column(column)
    column = table.Column(data=error,
                          name='eH77a_'+parname,
                          unit=unit)
    tbl.add_column(column)
 
# sort such that e10 comes after e9
import natsort
tbl = tbl[natsort.index_natsorted(tbl['ObjectName'])]

detection_note = ['-' if name in detections else
                  'weak' if name in weakdetections else
                  'none'
                  for name in tbl['ObjectName']]
tbl.add_column(table.Column(data=detection_note, name='DetectionStatus'))

ok = np.array([row['ObjectName'] in detections+weakdetections
               for row in tbl])

tbl[ok].write(paths.tpath('H77a_spectral_fits.ipac'), format='ascii.ipac')

for old,new in [('ObjectName','Object Name'),
                ('H77a_amplitude', 'Amplitude',),
                ('eH77a_amplitude','$\sigma$(Amplitude)',),
                ('H77a_velocity',  '$V_{LSR}$',),
                ('eH77a_velocity', '$\sigma(V_{LSR})$',),
                ('H77a_width',     '$dV (\sigma)$',),
                ('eH77a_width',    '$\sigma(dV)$'),
                ('DetectionStatus', 'Detection Status'),
               ]:
    tbl.rename_column(old, new)

for row in tbl:
    if "_" in row['Object Name']:
        row['Object Name'] = row['Object Name'].replace("_","\_")
latexdict['header_start'] = '\label{tbl:h77a}'
latexdict['caption'] = 'H$77\\alpha$ emission line parameters'
tbl[ok].write(paths.tpath('H77a_spectral_fits.tex'), format='ascii.latex', latexdict=latexdict)
