import paths
import ds9
from astropy.io import fits
from photom_files import files
import time

dd = ds9.ds9(target='w51points', start=True)
time.sleep(2)

names = ['25.0 GHz Epoch 2',
         '27.0 GHz Epoch 2',
         '2.5 GHz Epoch 2',
         '3.5 GHz Epoch 2',
         '4.9 GHz Epoch 3',
         '5.9 GHz Epoch 3',
         '8.4 GHz Epoch 1',
         '14.1 GHz Epoch 2',
         '12.6 GHz Epoch 2',
         '22.5 GHz Epoch 1',
         ]

dd.set('file {0}'.format(paths.dpath('Cband_Epoch3sm-Epoch3.fits')))
dd.set('frame new')
dd.set('file {0}'.format(paths.dpath('Kuband_Epoch3sm-Epoch3.fits')))

for name in names:
    dd.set('frame new')
    
    dd.set('file {0}'.format(files[name]))

    dd.set('wcs append', "OBJECT  = '{0}'".format(name))

dd.set('tile yes')
dd.set('frame frameno 1')
dd.set('wcs fk5')
dd.set('lock frame wcs')
