import numpy as np
import pylab as pl
from astropy.io import ascii
from astropy import table
from astropy.table import Table
import paths
from common_constants import distance
from astropy import units as u


tbl = Table.read(paths.tpath('EVLA_VLA_PointSourcePhotometry.ipac'), format='ascii.ipac')
lrad = (4*np.pi*distance**2 * tbl['peak_flux'] * u.beam * tbl['Frequency']).to(u.erg/u.s)

bieging = ascii.read(paths.tpath('bieging1989.tbl'))
beck = Table.read(paths.tpath('debecker2013.tbl'), format='ascii.commented_header')


pl.figure(1).clf()
pl.plot(bieging['mass'], np.transpose([bieging['f6cm_lo']*bieging['distance']**2,
                                       bieging['f6cm_hi']*bieging['distance']**2]),
        marker='o', linestyle='none')
pl.ylabel("$S_\\nu D^2$")
pl.xlabel("$M_\odot$")

pl.figure(2).clf()
pl.semilogx(10**bieging['luminosity'],
            np.transpose([bieging['f6cm_lo']*bieging['distance']**2,
                          bieging['f6cm_hi']*bieging['distance']**2]),
            marker='o', linestyle='none')
pl.ylabel("$S_\\nu D^2$")
pl.xlabel("$L_\odot$")

pl.figure(3).clf()
bins = np.linspace(-1.0,3.0,25)
fluxoverdist = np.log10(tbl['peak_flux'][tbl['Frequency']==4.9]*1000.*5.4**2)
pl.hist(fluxoverdist[np.isfinite(fluxoverdist)],
        bins=bins, histtype='stepfilled')
pl.hist(np.log10(bieging['f6cm_hi']*bieging['distance']**2), bins=bins,
        histtype='stepfilled', alpha=0.5)
pl.xlabel("$S_\\nu D^2$")

pl.figure(4).clf()
bins = np.linspace(27, 32)
pl.hist(np.log10(lrad.value), bins=bins, label='W51',
        histtype='step')
axis = pl.axis()
pl.hist(np.log10(beck['L_rad']), bins=bins,
        histtype='stepfilled',
        label='de Becker 2013', alpha=0.5)
pl.legend(loc='best')
pl.axis(axis)

pl.draw(); pl.show()
