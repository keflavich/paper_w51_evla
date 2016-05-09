"""
Results:
Mean velocity including only solid detections: 57.73833333333334
Velocity dispersion including only solid detections: 1.990338469261502
Virial mass in e1e2 using R=0.16 pc from spanning_tree: 1498.3488783876544 solMass
Density from that mass: 10^6.104362280824501 = 1.2716344379920037 1 / cm3x10^6
Mean velocity including weak detections: 58.05375
Velocity dispersion including weak detections: 4.298784530248055
Virial mass in e1e2 using R=0.16 pc (eyeballed): 6989.569498519303 solMass
Density from that mass: 10^6.773199761246572 = 5.93198113554171 1 / cm3x10^6
Mean velocity in e1: 58.086
Velocity dispersion in e1: 2.0071133500627205
Virial mass in e1 using R=0.07 pc (eyeballed): 666.6239757158603 solMass
Density from that mass: 10^6.829696091049973 = 6.756100341002692 1 / cm3x10^6
"""
import numpy as np
import paths
from astropy import table
from astropy import units as u
from astropy import constants


def mvir(radius, vdisp_1d, prefactor_eta=3.39):
    """
    Compute the virial mass given the 1-D velocity dispersion and the half-mass radius

    This is a simplified virial mass with
    M = sigma_3d**2 * radius / G
    and sigma_3d = sqrt(3) * sigma_1d

    prefactor_eta = 6*16/(3*pi)/3 = 3.39 for a Plummer sphere for stars
    prefactor_eta = 1 for r^-2 cloud (MacLaren 1988)
    prefactor_eta = 5/3 for rho=const cloud (MacLaren 1988)
    """
    return prefactor_eta * ((3 * vdisp_1d**2) * radius / constants.G).to(u.M_sun)

tbl22 = table.Table.read(paths.tpath('H2CO22_emission_spectral_fits.ecsv'), format='ascii.ecsv')
tbl77 = table.Table.read(paths.tpath('H77a_spectral_fits.ipac'), format='ascii.ipac')

cluster_sources = ['e8mol', 'e10mol', 'e1', 'e2', 'e3', 'e4', 'e8']
cluster_sources_optimistic = ['e8mol', 'e10mol', 'e1', 'e2', 'e3', 'e4', 'e8', 'e9', 'e10']
e1cluster = ['e8mol', 'e10mol', 'e1', 'e3', 'e4', 'e8',]

velocities = ([row['$V_{LSR}$'] for row in tbl22 if row['Object Name'] in cluster_sources] +
              [row['H77a_velocity'] for row in tbl77 if row['ObjectName'] in cluster_sources])
print("Mean velocity including only solid detections: {0}".format(np.mean(velocities)))
print("Velocity dispersion including only solid detections: {0}".format(np.std(velocities)))
print("Virial mass in e1e2 using R=0.16 pc from spanning_tree: {0}".format(mvir(0.16*u.pc, np.std(velocities)*u.km/u.s)))
density_e1e2 = (mvir(0.16*u.pc, np.std(velocities)*u.km/u.s)/(2.8*u.Da)/(4/3.*np.pi*(0.16*u.pc)**3)).to(u.cm**-3)
print("Density from that mass: 10^{0} = {1}x10^6".format(np.log10(density_e1e2.value), density_e1e2/1e6))
velocities_optimistic = ([row['$V_{LSR}$'] for row in tbl22 if row['Object Name'] in cluster_sources_optimistic] +
                         [row['H77a_velocity'] for row in tbl77 if row['ObjectName'] in cluster_sources_optimistic])
print()
print("Mean velocity including weak detections: {0}".format(np.mean(velocities_optimistic)))
print("Velocity dispersion including weak detections: {0}".format(np.std(velocities_optimistic)))
# spanning tree gives 0.09 pc
print("Virial mass in e1e2 using R=0.16 pc (eyeballed): {0}".format(mvir(0.16*u.pc, np.std(velocities_optimistic)*u.km/u.s)))
density_e1e2_wk = (mvir(0.16*u.pc, np.std(velocities_optimistic)*u.km/u.s)/(2.8*u.Da)/(4/3.*np.pi*(0.16*u.pc)**3)).to(u.cm**-3)
print("Density from that mass: 10^{0} = {1}x10^6".format(np.log10(density_e1e2_wk.value), density_e1e2_wk/1e6))
velocities_optimistic = ([row['$V_{LSR}$'] for row in tbl22 if row['Object Name'] in e1cluster] +
                         [row['H77a_velocity'] for row in tbl77 if row['ObjectName'] in e1cluster])

print()
print("Mean velocity in e1: {0}".format(np.mean(velocities_optimistic)))
print("Velocity dispersion in e1: {0}".format(np.std(velocities_optimistic)))
# spanning tree gives 0.04 pc
mvir_e1 = mvir(0.07*u.pc, np.std(velocities_optimistic)*u.km/u.s)
print("Virial mass in e1 using R=0.07 pc (eyeballed): {0}".format(mvir_e1))
density_e1 = (mvir_e1/(2.8*u.Da)/(4/3.*np.pi*(0.07*u.pc)**3)).to(u.cm**-3)
print("Density from that mass: 10^{0} = {1}x10^6".format(np.log10(density_e1.value), density_e1/1e6))

mvir_e1_gas = mvir(0.07*u.pc, np.std(velocities_optimistic)*u.km/u.s, prefactor_eta=5/3.)
print("Virial mass in e1 using R=0.07 pc (eyeballed) and gas eqn: {0}".format(mvir_e1_gas))
density_e1_gas = (mvir_e1_gas/(2.8*u.Da)/(4/3.*np.pi*(0.07*u.pc)**3)).to(u.cm**-3)
print("Density from that mass: 10^{0} = {1}x10^6".format(np.log10(density_e1_gas.value), density_e1_gas/1e6))

# Are the clusters star or gas dominated?
# see cluster_properties.py
# Implied stellar mass of e1 cluster, assuming 7 stars are OB stars: 695.874032382 +/- 238.642840322 M_sun
# Implied stellar mass of e1 cluster, assuming 7 stars are O stars: 2388.85335601 +/- 752.01450386 M_sun
