import numpy as np
from paths import dpath
from spectral_cube import SpectralCube
from astropy import units as u
from fnames import cube_names

errors = {}

for name, fn in cube_names.items():
    cube = SpectralCube.read(dpath(fn))
    print(cube)
    err = cube.std().to(u.mJy)
    print("{0} std: {1}".format(name, err))
    errors[name] = err
