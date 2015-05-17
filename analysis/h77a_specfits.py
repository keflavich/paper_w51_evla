""" was vdisp.py """
import pyspeckit 
import glob
from astropy import table
from astropy import units as u
from astropy import log
from astropy.utils.console import ProgressBar
import paths

sp = [pyspeckit.Spectrum(x) for x in
      ProgressBar(
          glob.glob(
              paths.dpath(
                  "spectra_h77/H77a_BDarray_speccube_uniform_contsub_cvel_big*fits")))]
spectra = sp

# conversion....
[s.xarr.convert_to_unit('km/s') for s in sp]

# setup
for s in sp:
    s.specname = s.fileprefix.split("_")[-1]
    log.info(s.specname+" stats")
    s.error[:] = s.stats((100,120))['std']

# fitting
for s in sp:
    log.info(s.specname+" fitting")
    s.plotter(xmin=-10,xmax=120)
    s.specfit(fittype='gaussian',
              guesses=[1,55,10],
              limited=[(True,False),(False,False),(True,False)])
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
    data = [sp.specfit.parinfo[ii].value
            for sp in spectra]
    error = [sp.specfit.parinfo[ii].error
            for sp in spectra]
    column = table.Column(data=data,
                          name='H77a_'+parname,
                          unit=unit)
    tbl.add_column(column)
    column = table.Column(data=error,
                          name='eH77a_'+parname,
                          unit=unit)
    tbl.add_column(column)

tbl.write(paths.tpath('H77a_spectral_fits.ipac'), format='ascii.ipac')
tbl.write(paths.tpath('H77a_spectral_fits.tex'), format='ascii.tex')
