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
import pylab as pl
pl.matplotlib.rc_file('pubfiguresrc')

h77a_freq = pyspeckit.spectrum.models.hydrogen.rrl(77)*u.GHz
he77a_freq = pyspeckit.spectrum.models.hydrogen.rrl(77, amu=4)*u.GHz
c77a_freq = pyspeckit.spectrum.models.hydrogen.rrl(77, amu=12)*u.GHz

# Read in all spectra extracted to the h77a directory matching the "best" h77a file
# May 18, 2015: I think H77a big2 is best, but I haven't inspected.  EDIT: big2 is in velocity.
# September 9, 2015: I think I deleted big2
# Exclude "mol" files since they are molecular and should be extracted from the molecular cubes
# Exclude "?" files since they are misnamed and use the obsolete scheme from
# before I decided upon e9/e10
sp = [pyspeckit.Spectrum(x) for x in
      ProgressBar(
          glob.glob(
              paths.dpath(
                  "spectra_h77/H77a_BDarray_speccube_uniform_contsub_cvel_big*e[0-9]*.fits")))
      if 'mol' not in x and '?' not in x
     ]
spectra = sp
# Also fit, but don't include in the table, the IRS 1 / IRS 2 / Lacy Jet spectra
# (important for inspecting for CRRLs)
sp_other = [pyspeckit.Spectrum(x) for x in
            ProgressBar(
                glob.glob(
                    paths.dpath(
                        "spectra_h77/H77a_BDarray_speccube_uniform_contsub_cvel_big_[ilr]*.fits")))
            if 'mol' not in x and '?' not in x
           ]

# My manual inspection: which are detected?
# weakdetections are those that are not clearly believable
detections = ['e1', 'e2', 'e3', 'e4', 'e6',]
weakdetections = ['e5', 'e9', 'e10']

# conversion....
[s.xarr.convert_to_unit('km/s') for s in sp+sp_other]

# setup
for ss in sp+sp_other:
    ss.data *= 1e3 # convery Jy->mJy
    ss.unit = 'mJy/beam'
    ss.specname = ss.fileprefix.split("_")[-1]
    log.info(ss.specname+" stats")
    noiseregion = (ss.xarr < 20*u.km/u.s).value | (ss.xarr > 100*u.km/u.s).value
    ss.error[:] = ss.data[noiseregion].std()

# fitting
for ii,ss in enumerate(sp+sp_other):
    fig = pl.figure(1)
    fig.clf()
    assert ss.xarr.unit == u.km/u.s
    ss.plotter(xmin=-10,xmax=120, figure=fig)
    try:
        assert ss.plotter.xlabel == 'Velocity (km / s)'
    except AssertionError:
        print("Label is {0} instead of Velocity (km / s)".format(ss.plotter.xlabel))
        print("Drawn label is {0}".format(ss.plotter.axis.get_xlabel()))
    ss.specfit(fittype='gaussian',
               guesses=[1,55,10],
               limited=[(True,False),(False,False),(True,False)])
    log.info(s.specname+" fitting: {0}".format(ss.specfit.parinfo))
    ss.plotter.ymin -= 0.3
    ss.plotter.label(ylabel='$S_\\nu$ (mJy beam$^{-1}$)')

    ax2 = ss.plotter.axis.twinx()
    ax2.set_ylim(*(np.array(ss.plotter.axis.get_ylim()) * ss.header['JYTOK']/1e3))
    ax2.set_ylabel("$T_B$ (K)")

    ss.specfit.plotresiduals(axis=ss.plotter.axis,clear=False,yoffset=-0.3,label=False)
    ss.plotter.savefig(paths.fpath('spectra/h77/'+ss.specname+"_h77a_fit.png"),
                                  bbox_inches='tight')

    pl.figure(ii).clf()
    ss.plotter(xmin=-10,xmax=120, figure=pl.figure(ii), axis=pl.subplot(2,1,1))
    ss.plotter.ymin -= 0.3
    ss.plotter.label(ylabel='$S_\\nu$ (mJy beam$^{-1}$)')
    ss.specfit.plot_fit()

    ax2 = ss.plotter.axis.twinx()
    ax2.set_ylim(*(np.array(ss.plotter.axis.get_ylim()) * ss.header['JYTOK']/1e3))
    ax2.set_ylabel("$T_B$ (K)")

    ss.specfit.plotresiduals(axis=ss.plotter.axis,clear=False,yoffset=-0.3,label=False)

    ss_he = ss.copy()
    ss_he.xarr.convert_to_unit(u.GHz)
    ss_he.xarr.refX = he77a_freq
    ss_he.xarr.convert_to_unit(u.km/u.s)
    ss_he.plotter(xmin=-10, xmax=120, axis=pl.subplot(2,1,2))
    ss_c = ss.copy()
    ss_c.xarr.convert_to_unit(u.GHz)
    ss_c.xarr.refX = c77a_freq
    ss_c.xarr.convert_to_unit(u.km/u.s)
    ss_c.plotter(xmin=-10, xmax=120, axis=pl.subplot(2,1,2), clear=False, color='blue')
    ss.plotter.savefig(paths.fpath('spectra/h77/'+ss.specname+"_h_he_c.png"),
                                  bbox_inches='tight')


tbl = table.Table()
names = table.Column(data=[sp.specname for sp in spectra], name='ObjectName')
tbl.add_column(names)
for ii,(parname,unit) in enumerate([('amplitude',u.mJy),
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
                ('eH77a_amplitude','$E$(Amplitude)',),
                ('H77a_velocity',  '$V_{LSR}$',),
                ('eH77a_velocity', '$E(V_{LSR})$',),
                ('H77a_width',     '$\sigma_V$',),
                ('eH77a_width',    '$E(\sigma_V)$'),
                ('DetectionStatus', 'Detection Status'),
               ]:
    tbl.rename_column(old, new)

for row in tbl:
    if "_" in row['Object Name']:
        row['Object Name'] = row['Object Name'].replace("_","\_")
latexdict['header_start'] = '\label{tab:h77a}'
latexdict['caption'] = 'H$77\\alpha$ emission line parameters'
tbl[ok].write(paths.tpath('H77a_spectral_fits.tex'), format='ascii.latex', latexdict=latexdict)
