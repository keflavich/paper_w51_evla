from astropy import units as u
from astropy.io import ascii

colnames = ['Object Name',
            'Amplitude',
            '$\sigma(Amplitude)$',
            '$V_{LSR}$',
            '$\sigma(V_{LSR})$',
            '$dV$',
            '$\sigma(dV)',
            '$\Omega_{ap}$',]
units = {'Amplitude':u.mJy.to_string(u.format.Latex),
         '$\sigma(Amplitude)$':u.mJy.to_string(u.format.Latex),
         '$V_{LSR}$':(u.km/u.s).to_string(u.format.Latex),
         '$\sigma(V_{LSR})$':(u.km/u.s).to_string(u.format.Latex),
         '$dV$':(u.km/u.s).to_string(u.format.Latex),
         '$\sigma(dV)':(u.km/u.s).to_string(u.format.Latex),
         '$\Omega_{ap}$':u.sr.to_string(u.format.Latex),}
latexdict = ascii.latex.latexdicts['AA']
latexdict['tabletype'] = 'table*'
latexdict['units'] = units
