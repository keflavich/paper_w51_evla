from astropy import units as u
from spectral_cube import SpectralCube
from paths import dpath
from fnames import cube_names
import ds9

dd = ds9.ds9()

cubes = {name: SpectralCube(dpath(fn)) for name,fn in cube_names.items()}

vcubes = {name: cube.with_spectral_unit(u.km/u.s) for name,cube in cubes.items()}

for name,vc in vcubes.items():
    vc.to_ds9(dd.....)
