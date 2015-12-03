from gaussfitter import gaussfit
from astropy import units as u
from spectral_cube import SpectralCube
import pyregion
import paths
from astropy.table import Table, Column
from astropy import wcs
import pylab as pl

regions = pyregion.open(paths.rpath('ch3oh_maser_spots_channellabels.reg'))

cube = SpectralCube.read(paths.dpath('ch3oh_256_e2zoom_chan550to700.image.fits'))
vcube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio',
                                rest_value=6.668518*u.GHz)

dx, dy = 4,4

parnames = ("height", "amplitude", "x", "y", "width_x", "width_y", "rota")
tbl = Table(names=['velocity'] + [x for p in parnames for x in p,'e'+p])

for ii in pl.get_fignums(): pl.figure(ii).clf()

for ii,region in enumerate(regions):
    channel = int(region.attr[1]['text'])

    x,y = cube.wcs.sub([wcs.WCSSUB_CELESTIAL]).wcs_world2pix([region.coord_list], 0)[0]

    # not sure why, but it looks like all the images are offset by 1 pix
    x = x+1
    y = y+1

    cutout = cube[channel, y-dy:y+dy, x-dx:x+dx].value

    rms = cube[channel, :20, :20].std()

    #params: (height, amplitude, x, y, width_x, width_y, rota)
    (pars, errors), fitimage = gaussfit(cutout, err=rms,
                                        params=(0, cutout[dy,dx], dx, dy, 1, 1, 0),
                                        limitedmin=[False,False,True,True,True,True,True],
                                        limitedmax=[False,False,True,True,False,False,True],
                                        minpars=[0,0,dx-1,dy-1,0,0,0],
                                        maxpars=[0,0,dx+1,dy+1,2.5,2.5,180],
                                        returnfitimage=True,
                                        return_error=True)
    #print(dict(zip(parnames, zip(pars, errors))))

    row = dict(zip(parnames, pars))
    row.update(dict(zip(['e'+p for p in parnames], errors)))
    row.update({'velocity': vcube.spectral_axis[channel]})
    tbl.add_row(row)
    tbl.pprint()

    pl.figure(ii / 16)
    ax = pl.subplot(4, 4, ii % 16 + 1)
    ax.imshow(cutout, cmap=pl.cm.gray_r)
    ax.contour(fitimage)
    ax.set_title(str(channel))

    #if ii / 16 > 0:
    #    break

tbl.write('CH3OH_fits.ipac', format='ascii.ipac')
pl.draw()
pl.show()
