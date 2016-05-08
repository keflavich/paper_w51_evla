from astroquery.vizier import Vizier
from astropy.table import Table

Vizier.ROW_LIMIT = 1e5
moxc = Vizier.query_constraints(catalog='J/ApJS/213/1/table3', GLON='>49.0 & < 51.0',)
w51moxc = moxc[0][(moxc[0]['Region'] == b'W51A')]
for cn in w51moxc.colnames:
    new = cn.replace("-","_")
    if cn != new:
        w51moxc.rename_column(cn, new)

w51moxc.write('../tables/w51_moxc.ecsv', format='ascii.ecsv')
w51moxc.write('../tables/w51_moxc.csv', format='ascii.csv')
w51moxc.write('../tables/w51_moxc.ipac', format='ascii.ipac')

# verify reading in what was written
Table.read('../tables/w51_moxc.csv', format='ascii.csv')
Table.read('../tables/w51_moxc.ipac', format='ascii.ipac')
#Table.read('../tables/w51_moxc.ecsv', format='ascii.ecsv')

with open('../regions/moxc.reg', 'w') as f:
    f.write("fk5\n")
    for row in w51moxc:
        f.write("point({0},{1}) # point=x color=red text={{{2}}}\n".format(row['RAJ2000'], row['DEJ2000'], row['CXOU']))


import numpy as np
import pylab as pl
from wcsaxes import WCS as WCSaxes
from astropy import wcs
from astropy.io import fits
import paths
from astropy.convolution import convolve, Gaussian2DKernel

hdu = fits.open(paths.dpath("W51C_ACarray_continuum_4096_both_uniform_contsplit.clean.image.fits"))[0]
mywcs = wcs.WCS(hdu.header).sub([wcs.WCSSUB_CELESTIAL])
wcsaxes = WCSaxes(mywcs.to_header())

fig = pl.figure(1)
fig.clf()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcsaxes)
im = ax.imshow(hdu.data.squeeze()*1e3, cmap=pl.cm.gray_r, origin='lower')
ra = ax.coords['ra']
ra.set_major_formatter('hh:mm:ss.s')
dec = ax.coords['dec']
ra.set_axislabel("RA (J2000)")
dec.set_axislabel("Dec (J2000)")
lims = ax.axis()
mylims = ((290.94465,  14.492588), (290.90729,14.530606))
(x1,y1),(x2,y2) = mywcs.wcs_world2pix(mylims, 0)

clims = mywcs.wcs_pix2world([[lims[0],lims[2]], [lims[1],lims[3]]], 0)
bins = [np.linspace(clims[1,0], clims[0,0], 200),
        np.linspace(clims[0,1], clims[1,1], 200)]

tr_fk5 = ax.get_transform("fk5")

dots, = ax.plot(w51moxc['RAJ2000'], w51moxc['DEJ2000'], 'r.', transform=tr_fk5,
                markersize=2)
#ax.axis(lims)
ax.axis([x1,x2,y1,y2])

fig.savefig(paths.fpath("moxc_points_on_cband.png"))

H,bx,by = np.histogram2d(w51moxc['RAJ2000'], w51moxc['DEJ2000'], bins=bins)
H2 = convolve(H, Gaussian2DKernel(2))
cx = (bx[1:]+bx[:-1])/2.
cy = (by[1:]+by[:-1])/2.
con = ax.contour(cx,
                 cy,
                 H2.T, transform=tr_fk5,
                 levels=[0.12,  0.18,  0.24,  0.3 ,  0.36,  0.42],
                 colors=['b']*10,
                 linewidth=0.5,
                 interpolation='bicubic',
                )
#ax.axis(lims)
ax.axis([x1,x2,y1,y2])
fig.savefig(paths.fpath("moxc_points_contours_on_cband.png"))
dots.set_visible(False)
fig.savefig(paths.fpath("moxc_contours_on_cband.png"))
pl.draw()
pl.show()

