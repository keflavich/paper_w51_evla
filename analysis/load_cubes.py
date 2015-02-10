from astropy import units as u
from spectral_cube import SpectralCube
from paths import dpath
from fnames import cube_names, continua
import ds9

dd = ds9.ds9()

cubes = {name: SpectralCube.read(dpath(fn)) for name,fn in cube_names.items()}

vcubes = {name: cube.with_spectral_unit(u.km/u.s, velocity_convention='radio')
          for name,cube in cubes.items()}

for name,vc in vcubes.items():
    vc.header['OBJECT'] = name
    vc.to_ds9(dd.id, newframe=True)

dd.set('tile yes')
dd.set('frame frameno 1')
dd.set('frame delete')
dd.set('scale limits -0.005 0.005')
dd.set('wcs fk5')
dd.set('lock frame wcs')
dd.set('lock slice wcs')
dd.set('lock scale yes')
dd.set('lock color yes')

for name,cfn in continua.items():
    dd.set('frame new')
    f = fits.open(dpath(cfn))
    f[0].header['OBJECT'] = name
    dd.set_pyfits(f)
