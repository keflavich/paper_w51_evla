""" never tested, never worked. Probably could work if I got the right keywords... """
from astropy import units as u
from astropy import coordinates
from astroquery.simbad import Simbad

irs2 = coordinates.SkyCoord.from_name('W51 IRS2')

region = 'region(circle, fk5, {ra} {dec}, {rad_deg})'.format(ra=irs2.ra.to(u.deg).value,
                                                             dec=irs2.dec.to(u.deg).value,
                                                             rad_deg=12./60.)

result = Simbad.query_criteria(region, otype='star', sptype='OB')
