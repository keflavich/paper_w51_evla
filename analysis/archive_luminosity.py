from astropy import units as u
from astropy import coordinates
from astroquery.irsa import Irsa
from common_constants import distance

IRAS = Irsa.query_region(coordinates.SkyCoord.from_name('W51'),
                         catalog='iraspsc', radius=1*u.arcmin)
Akari = Irsa.query_region(coordinates.SkyCoord.from_name('W51'),
                          catalog='akari_fis', radius=1*u.arcmin)

# formulae from http://marc.sauvage.free.fr/astro_book/IRAS_pages/IRAS.html
fir_lum_iras = 3.96e5 * (2.58*IRAS['fnu_60'][0] + IRAS['fnu_100'][0]) * (distance.to(u.Mpc).value)**2 * u.L_sun
mir_lum_iras = 1.611e6 * (2.61*IRAS['fnu_12'][0] + IRAS['fnu_25'][0]) * (distance.to(u.Mpc).value)**2 * u.L_sun

print("IRAS FIR luminosity: {0}".format(fir_lum_iras))
print("IRAS MIR luminosity: {0}".format(mir_lum_iras))
print("Harvey 1986 luminosity: {0}".format(1e7*u.L_sun))
print("Sievers 1991 luminosity: {0}".format(1.8e7*u.L_sun*(7.5*u.kpc/distance)**-2))

# Harvey 1986a 1986ApJ...300..737H: 10^7 Lsun

import pylab as pl

iras_wl = [12,25,60,100]
akari_wl = [65,90,140,160]


pl.figure(1).clf()
pl.plot(iras_wl, [IRAS['fnu_{0}'.format(wl)] for wl in iras_wl], 's')
pl.plot(akari_wl, [Akari['flux{0}'.format(wl)] for wl in akari_wl], 'o')
