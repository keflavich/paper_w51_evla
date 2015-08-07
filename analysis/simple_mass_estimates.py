from __future__ import print_function
import numpy as np
import pyradex
from astropy import units as u
from paths import tpath
import paths
from astropy.table import Table, Column

R = pyradex.Radex(species='oh2co-h2', density=1e6, abundance=1e-9, temperature=30*u.K)
print(R(density=10**6.5, abundance=1e-9, temperature=35*u.K)['T_B','tau','Tex'][2])

abundance=1e-9
# Compute brightness as a function of density for fixed abundance for two temperatures
pars = ('T_B','tau','Tex')
if 'results' not in locals():
    results = {tem:
               {dens:
                {k:v for k,v in zip(pars, R(abundance=abundance,
                                            collider_densities={'H2':dens},
                                            temperature=tem)[pars][2])}
                for dens in np.logspace(3,8,100)
               }
               for tem in [30,60]
              }

emi_tbl = Table.read(paths.tpath('H2CO22_emission_spectral_fits.ecsv'), format='ascii.ecsv')
measured_brightness = emi_tbl['Amplitude'] * 32.398

tb_60 = np.array([row['T_B'] for row in results[60].values()])
bestdens_60 = [results[60].keys()[np.argmin(np.abs(tb_60-mb))]
               for mb in measured_brightness]*u.cm**-3
print("Density estimates for T=60: ",np.log10(bestdens_60.value))
tb30 = np.array([row['T_B'] for row in results[30].values()])
bestdens30 = [results[30].keys()[np.argmin(np.abs(tb30-mb))]
               for mb in measured_brightness]*u.cm**-3
print("Density estimates for T=30: ",np.log10(bestdens30.value))


distance = 5.41*u.kpc
linewidth = 5*u.km/u.s
sourcearea = (emi_tbl['$\Omega_{ap}$'] * distance**2).to(u.cm**2)
sourcesize = sourcearea**0.5
particle_mass = 2.8*u.Da

source_mass = (((4/3.*np.pi*sourcesize**3) * particle_mass * (bestdens_60)).to(u.M_sun))
print(source_mass)

emi_tbl.add_column(Column(data=source_mass, name="$M_{75 K}$"))
emi_tbl.add_column(Column(data=bestdens_60, name="$n_{75 K}(H_2)"))

#source_mass_v2 = ((sourcearea * particle_mass * ).to(u.M_sun))
#print("{0:0.3g}".format(source_mass))
