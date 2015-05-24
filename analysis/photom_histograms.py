"""
Run ptsrc_photom interactively first
"""

import pylab as pl
import numpy as np
import astroML.plotting as ampl


pointsources = [
 'e1c?',
 'AA',
 'e1d?',
 'AC',
 'AB',
 'AE',
 'AD',
 'AG',
 'AF',
 'e1b?',
 'e5',
 'e4',
 'e6',
 'e1',
 'e3',
 'e2',
 ]

pl.figure(2)
pl.clf()
for ii,frq in enumerate(fluxes.keys()):
    pl.subplot(3,3,ii+1)
    d = np.array([fluxes[frq][x] for x in pointsources])
    ampl.hist(3+np.log10(d[d>0]), bins=10, log=True, alpha=0.5, histtype='stepfilled')
    pl.title(frq)

pl.figure(3)
pl.clf()
for ii,frq in enumerate(peaks.keys()):
    pl.subplot(3,3,ii+1)
    d = np.array([peaks[frq][x] for x in pointsources])
    ampl.hist(3+np.log10(d[d>0]), bins=10, log=True, alpha=0.5, histtype='stepfilled')
    pl.title(frq)

pl.figure(4)
pl.clf()
frq = '14.1 GHz'

d = np.array([peaks[frq][x] for x in pointsources])
ampl.hist(3+np.log10(d[d>0]), bins=np.linspace(-0.55,2.01,5), log=True, alpha=0.5, histtype='stepfilled')
pl.title(frq)

pl.show()
