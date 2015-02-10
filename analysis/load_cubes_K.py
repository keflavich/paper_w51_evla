from astropy import units as u
from spectral_cube import SpectralCube
from paths import dpath
from fnames import cube_names, continua
import radio_beam
import ds9

dd = ds9.ds9()

cubes = {name: SpectralCube.read(dpath(fn)) for name,fn in cube_names.items()}

vcubes = {name: cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
          for name,cube in cubes.items()}

beams = {name: radio_beam.Beam.from_fits_header(cube.header)
         for name,cube in vcubes.items()}

for cube in vcubes.values():
    cube._unit = u.Jy

vkcubes = {name: cube.with_flux_unit(u.K,
                                     u.brightness_temperature(beams[name],
                                                              cube.wcs.wcs.restfrq*u.Hz))
          for name,cube in vcubes.items()}

for name,vc in vkcubes.items():
    vc._header['OBJECT'] = name
    vc.to_ds9(dd.id, newframe=True)

dd.set('tile yes')
dd.set('frame frameno 1')
dd.set('frame delete')
dd.set('scale limits -300 300')
dd.set('wcs fk5')
dd.set('lock frame wcs')
dd.set('lock slice wcs')
dd.set('lock scale yes')
dd.set('lock color yes')

for name,cfn in continua.items():
    dd.set('frame new')
    f = fits.open(dpath(cfn))
    beam = radio_beam.Beam.from_fits_header(f[0].header)
    freq = 4.82966*u.GHz if '11' in name else 14.488*u.GHz
    JyToK = (1*u.Jy).to(u.K, u.brightness_temperature(beam, freq)).value
    f[0].data = f[0].data*JyToK
    f[0].header['OBJECT'] = name
    dd.set_pyfits(f)
