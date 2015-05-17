import numpy as np
import pyspeckit
from astropy import units as u
import radio_beam

bm = radio_beam.Beam.from_fits_header('W51-36GHzcont.map.image.fits')
x,y=np.loadtxt('spec-w51e2-tot-66-pb.txt', comments="#").T
amf = pyspeckit.spectrum.models.ammonia.ammonia_model_background()
xarr = pyspeckit.spectrum.units.SpectroscopicAxis(x*u.km/u.s,
                                                  refX=pyspeckit.spectrum.models.ammonia.freq_dict['sixsix']*u.Hz,
                                                  velocity_convention='radio')
sp = pyspeckit.Spectrum(xarr=xarr, data=y)
spK = sp.copy()
spK.data = sp.data * ((1*u.Jy).to(u.K, u.brightness_temperature(bm,
                                                                sp.xarr.refX)).value)
spK.specfit.Registry.add_fitter('ammonia_bg', amf, 7)

# Fit each component
sp.specfit(fittype='gaussian', guesses=[-0.15, 25, 3,
                                        -0.15, 32, 3,
                                        -0.30, 58, 6,
                                        -0.15, 84, 3,
                                        -0.15, 90, 3])
print "Velocity Offsets: ",[(((x-sp.specfit.parinfo[7].value)*u.km/u.s)) for x in sp.specfit.parinfo.values[1::3]]
print "Frequency Offsets: ",[((((x-sp.specfit.parinfo[7].value)*u.km/u.s)/constants.c) * sp.xarr.refX).to(u.MHz) for x in sp.specfit.parinfo.values[1::3]]
models = pyspeckit.spectrum.models.ammonia_constants.voff_lines_dict['sixsix'][:3] + pyspeckit.spectrum.models.ammonia_constants.voff_lines_dict['sixsix'][-2:]
print "Difference between observed and modeled offsets: ",[(((x-sp.specfit.parinfo[7].value-m)*u.km/u.s)) for x,m in zip(sp.specfit.parinfo.values[1::3], models[::-1])]

T = True
F = False
#spK.specfit(fittype='ammonia_bg',
#            guesses=[300, 20000, 20, 1.34, 58, 0.5, 4.2e4],
#            fixed=[T,F,F,F,F,T,T],
#            verbose=True, quiet=False, shh=False, maxiter=1)
parinfo = amf.make_parinfo(params=[300, 20000, 20, 1.34, 58, 0.5, 4.2e4],
                            fixed=[F,F,F,F,F,T,T],)
parinfo.background_tb0.limited = (False,False)
parinfo.tex0.mpmaxstep=1000
print parinfo

spK.plotter()
spK.specfit.selectregion(exclude=[40,70], highlight=True)
spK.specfit(fittype='ammonia_bg', parinfo=parinfo, verbose=True, quiet=False,
            shh=False, veryverbose=True, reset_selection=False)

spK.specfit.plot_fit()
spK.specfit.plotresiduals(axis=spK.plotter.axis, clear=False, yoffset=-22000, label=False)
spK.plotter.axis.set_ylim(-25000, 5000)
spK.plotter.savefig("nh3_66_hyperfineonly.png")

spK.plotter()
spK.specfit(fittype='ammonia_bg', parinfo=parinfo, verbose=True, quiet=False,
            shh=False, veryverbose=True, reset_selection=True)
spK.specfit.plotresiduals(axis=spK.plotter.axis, clear=False, yoffset=-22000, label=False)
spK.plotter.axis.set_ylim(-25000, 5000)
spK.plotter.savefig("nh3_66_fitwhole.png")

parinfo = amf.make_parinfo(params=[300, 20000, 20, 1.34, 58, 0.5, 4.2e4,
                                   300, 20000, 19.5, 0.5, 59, 0.5, 4.2e4],
                           fixed=[T,F,F,F,F,T,T]*2,
                           npeaks=2)
parinfo.background_tb0.limited = (False,False)
parinfo.background_tb1.limited = (False,False)
parinfo.tex0.mpmaxstep=1000
parinfo.tex1.mpmaxstep=1000
print parinfo
spK.plotter()
spK.specfit(fittype='ammonia_bg', parinfo=parinfo, verbose=True, quiet=False,
            shh=False, veryverbose=True, reset_selection=True)
spK.specfit.plotresiduals(axis=spK.plotter.axis, clear=False, yoffset=-22000, label=False)
spK.plotter.axis.set_ylim(-25000, 5000)
spK.plotter.axis.set_xlim(-30,200)
spK.plotter.savefig("nh3_66_fitwhole_twocomp.png")

import pylab as pl
pl.draw()
pl.show()
