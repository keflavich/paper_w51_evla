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
units = {'Amplitude':u.mJy.to_string(u.format.LatexInline),
         '$\sigma(Amplitude)$':u.mJy.to_string(u.format.LatexInline),
         '$V_{LSR}$':(u.km/u.s).to_string(u.format.LatexInline),
         '$\sigma(V_{LSR})$':(u.km/u.s).to_string(u.format.LatexInline),
         '$dV$':(u.km/u.s).to_string(u.format.LatexInline),
         '$\sigma(dV)':(u.km/u.s).to_string(u.format.LatexInline),
         '$\Omega_{ap}$':u.sr.to_string(u.format.LatexInline),}
latexdict = ascii.latex.latexdicts['AA']
latexdict['tabletype'] = 'table*'
latexdict['tablealign'] = 'htp'
latexdict['units'] = units
