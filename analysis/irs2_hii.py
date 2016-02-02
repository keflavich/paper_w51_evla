import numpy as np
import pyspeckit
from common_constants import distance
from astropy import units as u
from astropy import constants
from spectral_cube import SpectralCube
import paths

cube = SpectralCube.read(paths.dpath('H77a_BDarray_speccube_briggs0_contsub_cvel_big.fits'))
h77ajytok = cube.beam.jtok(cube.wcs.wcs.restfrq*u.Hz)

# approximately 10 mJy/beam
# hmm... or 100 mJy/beam/(km/s)
#peak_tb_irs2 = 3500*u.K/(u.km/u.s)
peak_tb_irs2 = 100*u.mJy/(u.km/u.s) * h77ajytok/u.Jy
Te = 1e4*u.K
dnu = (10*u.km/u.s / constants.c * pyspeckit.spectrum.models.hydrogen.rrl(77)*u.GHz).to(u.kHz)
r_irs2_hii = (1.5*u.arcsec*distance).to(u.pc, u.dimensionless_angles())

# eqn 14.28 of Wilson 2009
EM_IRS2 = ((peak_tb_irs2/(u.K/(u.km/u.s))) / 1.92e3 * ((Te/u.K)**1.5) * (dnu/u.kHz) * u.cm**-6 * u.pc).to(u.cm**-6*u.pc)

n_IRS2 = ((EM_IRS2 / (2*r_irs2_hii))**0.5).to(u.cm**-3)
M_IRS2 = (4/3.*np.pi*r_irs2_hii**3 * n_IRS2 * 1.4 * constants.m_p).to(u.M_sun)

print("IRS2 EM={EM:0.2g} n={dens:0.2g}  M={mass:0.2g}".format(EM=EM_IRS2, dens=n_IRS2, mass=M_IRS2))

peak_tb_irs2outflow = 1.7*u.mJy/(u.km/u.s) * h77ajytok/u.Jy
r_irs2outflow = (9.5*u.arcsec*distance).to(u.pc, u.dimensionless_angles())

EM_irs2outflow = ((peak_tb_irs2outflow/(u.K/(u.km/u.s))) / 1.92e3 * ((Te/u.K)**1.5) * (dnu/u.kHz) * u.cm**-6 * u.pc).to(u.cm**-6*u.pc)

n_irs2outflow = ((EM_irs2outflow / (2*r_irs2outflow))**0.5).to(u.cm**-3)
M_irs2outflow = (4/3.*np.pi*r_irs2outflow**3 * n_irs2outflow * 1.4 * constants.m_p).to(u.M_sun)

print("IRS2 outflow EM={EM:0.2g} n={dens:0.2g}  M={mass:0.2g}".format(EM=EM_irs2outflow, dens=n_irs2outflow, mass=M_irs2outflow))

# IRS2 outflow is going through shells at different velocities
# We look at the furthest shell at 37.5 km/s (from section 3.1)
# The IRS2 center/peak velocity is 62.5 km/s
shell_width = (1.0*u.arcsec*distance).to(u.pc, u.dimensionless_angles())
vshell = 37.5*u.km/u.s
virs2 = 62.6*u.km/u.s
d_irs2_shell = (7.5*u.arcsec*distance).to(u.pc, u.dimensionless_angles())

volume_shell = (4/3.*np.pi*(r_irs2outflow**3 - (r_irs2outflow-shell_width)**3))
n_shell = ((EM_irs2outflow / (shell_width))**0.5).to(u.cm**-3)
M_shell = (volume_shell * n_irs2outflow * 1.4 * constants.m_p).to(u.M_sun)

print("IRS2 outflow shell EM={EM:0.2g} n={dens:0.2g}  M={mass:0.2g}".format(EM=EM_irs2outflow, dens=n_shell, mass=M_shell))

vshell = np.abs(virs2-vshell)
shell_timescale = (d_irs2_shell / vshell).to(u.kyr)

# upper limit from assuming the circle at 37.5 km/s is completely filled
upper_limit_massloss_rate = M_irs2outflow / shell_timescale

print("Timescale for IRS2 shell to reach location at current velocity v={vel} is "
      "t={timescale}, mass loss rate total "
      "mdot={mlrate}".format(timescale=shell_timescale,
                             mlrate=upper_limit_massloss_rate,
                             vel=vshell))
mclump = 1e4*u.M_sun
tevap = (mclump/upper_limit_massloss_rate).to(u.Myr)
print("Evaporation timescale t={tevap}".format(tevap=tevap))


rclump = (1.5*u.pc)
clump_density = mclump / (4/3. * np.pi * rclump**3)
freefalltime = ((clump_density*constants.G)**-0.5).to(u.Myr)
print("Free fall time for an M={mclump:0.2g} cluster with density n={clump_dens:0.2g} is tff={tff:0.2g}"
      .format(tff=freefalltime, mclump=mclump,
              clump_dens=(clump_density/(2.4*constants.m_p)).to(u.cm**-3)))

print("evaporation / free fall = {nfreefalls:0.4g}".format(nfreefalls=(tevap/freefalltime).decompose()))

print("escape speed: v={vesc:0.3g}".format(vesc=((2*constants.G*mclump/rclump)**0.5).to(u.km/u.s)))
