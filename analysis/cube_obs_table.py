import numpy as np
import paths
from paths import dpath
from spectral_cube import SpectralCube
from astropy import units as u
from fnames import cube_names,rrl_cube_names
from astropy.table import Table,Column
from latex_info import latexdict, format_float

cube_names.update(rrl_cube_names)

errors = {}
cubes = {}
beams = {}
dv = {}

for name, fn in cube_names.items():
    cube = SpectralCube.read(dpath(fn)).with_spectral_unit(u.km/u.s,
                                                           velocity_convention='radio')
    print(cube)
    err = cube.std().to(u.mJy)
    print("{0} std: {1}".format(name, err))
    errors[name] = err
    cubes[name] = cube
    beams[name] = cube.beam
    dv[name] = np.diff(cube.spectral_axis).mean()

name_mapping = {'11_natural': '1-1 Natural', '11_uniform': '1-1 Uniform',
                '22_briggs0': '2-2 Briggs 0', '22_natural': '2-2 Natural',
                'H77a': 'H$77\\alpha$ Briggs 0'}

tbl = Table([Column(data=[name_mapping[name] for name in cube_names],
                    name="Cube ID"),
             Column(data=([cubes[name].wcs.wcs.restfrq for name in cube_names]*u.Hz).to(u.GHz),
                    name='Frequency'),
             Column(data=[dv[name].to(u.km/u.s).value for name in cube_names]*u.km/u.s,
                    name='Channel Width'),
             Column(data=[errors[name].to(u.mJy).value for name in
                          cube_names]*u.mJy/u.beam, name='RMS'),
             Column(data=[beams[name].major.to(u.arcsec).value for name in
                          cube_names]*u.arcsec, name='BMAJ'),
             Column(data=[beams[name].minor.to(u.arcsec).value for name in
                          cube_names]*u.arcsec, name='BMIN'),
             Column(data=[beams[name].pa.to(u.deg).value for name in
                          cube_names]*u.deg, name='BPA'),

            ]
           )
tbl.add_column(Column(data=[(1*u.Jy).to(u.K, u.brightness_temperature(bm.sr,
                                                                      freq*u.GHz,)).value
                            for bm,freq in zip(beams.values(),
                                               tbl['Frequency'])
                           ],
                      name='Jy-Kelvin'))
tbl.add_column(Column(data=((tbl['RMS']*u.beam).to(u.Jy).value*tbl['Jy-Kelvin'])*u.K,
                      name='RMS$_K$'), index=4)

tbl.sort(['Frequency'])
latexdict['header_start'] = '\label{tab:cubes}'
latexdict['caption'] = 'Spectral Cubes'
latexdict['tablefoot'] = ('\par\nJy-Kelvin gives the conversion factor from Jy '
                          'to Kelvin given the synthesized beam size and '
                          'observation frequency.')
#latexdict['col_align'] = 'lllrr'
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
def round_to_n(x, n):
    return round(x, -int(np.floor(np.log10(x))) + (n - 1))


tbl.write(paths.tpath('cubes.tex'), format='ascii.latex',
          latexdict=latexdict,
          formats={'BMAJ': lambda x: '{0:0.2f}'.format(x),
                   'BMIN': lambda x: '{0:0.2f}'.format(x),
                   'BPA':  lambda x: '{0:0.2f}'.format(x),
                   'RMS': lambda x: '{0:0.2f}'.format(x),
                   'RMS$_K$': format_float,
                   'Channel Width': lambda x: '{0:0.2f}'.format(x),
                   'Frequency': lambda x: '{0:0.6f}'.format(x),
                   'Jy-Kelvin':  format_float,
                  },)
