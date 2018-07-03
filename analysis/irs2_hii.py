import numpy as np
import pyspeckit
from common_constants import distance
from astropy import units as u
from astropy import constants
from astropy import wcs
from spectral_cube import SpectralCube
import paths
import dust_emissivity
from astropy.io import fits
import radio_beam
import regions
import imf

cube = SpectralCube.read(paths.dpath('H77a_BDarray_speccube_briggs0_contsub_cvel_big.fits'))
h77ajytok = cube.beam.jtok(cube.wcs.wcs.restfrq*u.Hz)

# approximately 10 mJy/beam
# hmm... or 100 mJy/beam/(km/s)
#peak_tb_irs2 = 3500*u.K/(u.km/u.s)
# I don't really know where this came from...
#peak_tb_irs2 = 100*u.mJy/(u.km/u.s) * h77ajytok/u.Jy
# slice 200 contains the overall peak value, and the peak is on IRS2
peak_tb_irs2 = u.Quantity(cube[200,:,:].max())/(u.km/u.s) * h77ajytok/u.Jy
# if I drop an aperture on IRS2 and take the mean, it peaks at ~0.007 Jy
mean_tb_irs2 = 7*u.mJy/(u.km/u.s) * h77ajytok/u.Jy
Te = 1e4*u.K
dnu = (10*u.km/u.s / constants.c * pyspeckit.spectrum.models.hydrogen.rrl(77)*u.GHz).to(u.kHz)
r_irs2_hii = (1.5*u.arcsec*distance).to(u.pc, u.dimensionless_angles())

# eqn 14.28 of Wilson 2009
EM_IRS2 = ((mean_tb_irs2/(u.K/(u.km/u.s))) / 1.92e3 * ((Te/u.K)**1.5) * (dnu/u.kHz) * u.cm**-6 * u.pc).to(u.cm**-6*u.pc)

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
M_shell = (volume_shell * n_shell * 1.4 * constants.m_p).to(u.M_sun)

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
print("4pi r^2 v * n_shell = {0}".format((1.4*u.Da*n_shell*4*np.pi*r_irs2outflow**2*vshell).to(u.Msun/u.yr)))
print("4pi r^2 v * n_outflow = {0}".format((1.4*u.Da*n_irs2outflow*4*np.pi*r_irs2outflow**2*vshell).to(u.Msun/u.yr)))
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

print("\nSome extra calculations that aren't really useful.")
EM_irs2 = ((peak_tb_irs2/(u.K/(u.km/u.s))) / 1.92e3 * ((Te/u.K)**1.5) * (dnu/u.kHz) * u.cm**-6 * u.pc).to(u.cm**-6*u.pc)
n_irs2 = ((EM_irs2 / (2*r_irs2_hii))**0.5).to(u.cm**-3)
print("Density of the IRS2 HII region: {0}".format(n_irs2))
freefalltime_irs2hii = ((n_irs2*1.37*constants.m_p*constants.G)**-0.5).to(u.Myr)
print("Free-fall time of something at IRS2 HII's density: {0}".format(freefalltime_irs2hii))


print()
print("Actual calculations based on the ~2000 Msun of gas that's definitely star-forming")
almacont = fits.open(paths.almadpath('w51_te_continuum_best.fits'))
almabm = radio_beam.Beam.from_fits_header(almacont[0].header)
almawcs = wcs.WCS(almacont[0].header)
irs2reg = regions.io.ds9.read.DS9Parser('fk5; ellipse(19:23:39.912,+14:31:05.381,3.352",1.254",3.5839317e-06)').shapes.to_regions()[0]
irs2preg = irs2reg.to_pixel(almawcs)
irs2mask = irs2preg.to_mask()
total_flux = (irs2mask.cutout(almacont[0].data) * irs2mask.data).sum() * u.Jy
ppbeam = (almabm.sr / (wcs.utils.proj_plane_pixel_area(almawcs)*u.deg**2)).decompose()
total_flux_Jy = total_flux / ppbeam
m40 = dust_emissivity.dust.massofsnu(snu=total_flux_Jy, temperature=40*u.K, nu=226*u.GHz, distance=distance)
print('total flux: {0}  mass(40K) = {1}'.format(total_flux_Jy, m40))
reff = (((irs2reg.width * irs2reg.height / 4 * np.pi))**0.5 * distance).to(u.pc, u.dimensionless_angles())
vol = 4/3 * np.pi * reff**3
dens = (m40/vol)
tff = ((3 * np.pi / (32 * constants.G * dens))**0.5).decompose()
print("tff = {0}  (dens={1})".format(tff.to(u.yr), (dens/(2.8*u.Da)).to(u.cm**-3)))

mass_lost_per_tff = (upper_limit_massloss_rate * tff).to(u.M_sun)
print("total mass lost during 1 tff = {0}".format(mass_lost_per_tff))
print("total mass lost as fraction of total mass = {0}".format(mass_lost_per_tff/m40))
print("total mass lost as fraction of total mass * eps_ff=0.01 = {0}".format(mass_lost_per_tff/(m40*0.01)))

# kim kim ostriker #'s
# email of June 4, 2018: dot N_e = c_i [ Q R/alpha_B]^1/2
c_i = 10*u.km/u.s
alpha_B = 2.6e-13*u.cm**3*u.s
N_e_dot = (upper_limit_massloss_rate / (2.8*u.Da)).to(u.s**-1)
Q = ((N_e_dot / c_i)**2 * alpha_B / r_irs2outflow).decompose()
print("Q_lyc needed for our inferred mass loss rate at the radius of the outflow: {0} (log={1})".format(Q, np.log10(Q.value)))
print("Typical Q_lyc(5000msun cluster) = {0}".format(np.mean([imf.lyc_of_cluster(imf.make_cluster(5000, silent=True)) for x in range(20)])))
print("Typical Q_lyc(2000msun cluster) = {0}".format(np.mean([imf.lyc_of_cluster(imf.make_cluster(2000, silent=True)) for x in range(20)])))
print("Typical Q_lyc(10000msun cluster) = {0}".format(np.mean([imf.lyc_of_cluster(imf.make_cluster(10000, silent=True)) for x in range(20)])))
