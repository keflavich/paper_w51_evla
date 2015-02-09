from astropy import units as u
from spectral_cube import SpectralCube
from paths import dpath
from fnames import cube_names
import ds9

dd = ds9.ds9()

cubes = {name: SpectralCube.read(dpath(fn)) for name,fn in cube_names.items()}

vcubes = {name: cube.with_spectral_unit(u.km/u.s) for name,cube in cubes.items()}

for name,vc in vcubes.items():
    vc.to_ds9(dd.id, newframe=True)

dd.set('tile yes')
dd.set('frame frameno 1')
dd.set('frame delete')
dd.set('scale limits -0.1 0.1')
dd.set('wcs fk5')
dd.set('lock frame wcs')
dd.set('lock slice wcs')
dd.set('lock scale yes')
dd.set('lock color yes')
