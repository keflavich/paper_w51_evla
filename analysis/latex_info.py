from astropy import units as u
from astropy.io import ascii

def exp_to_tex(st):
    if st == 'nan':
        return '-'
    elif 'e' in st:
        pt1,pt2 = st.split('e')
        return "{0}\\ee{{{1:d}}}".format(pt1,int(pt2))
    return st

def format_float(st):
    return exp_to_tex("{0:0.2g}".format(st))


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
