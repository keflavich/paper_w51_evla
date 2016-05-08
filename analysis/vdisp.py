"""
Results:
Mean velocity including only solid detections: 57.428999999999995
Velocity dispersion including only solid detections: 1.967670958265126
Virial mass in e1e2 using R=0.16 pc from spanning_tree: 431.9806822908472 solMass
Density from that mass: 366617.8952524318 1 / cm3
Mean velocity including weak detections: 57.87785714285714
Velocity dispersion including weak detections: 4.544453728874081
Virial mass in e1e2 using R=0.16 pc from spanning_tree: 2304.2136568113438 solMass
Density from that mass: 1955564.2085479042 1 / cm3
Mean velocity in e1: 57.78625
Velocity dispersion in e1: 2.0497618001855726
Virial mass in e1 using R=0.07 pc from spanning_tree: 205.08987600154953 solMass
Density from that mass: 2078544.7743644696 1 / cm3
"""
import numpy as np
import paths
from astropy import table
from astropy import units as u
from astropy import constants


def mvir(radius, vdisp):
    """
    Compute the virial mass given the 1-D velocity dispersion and the half-mass radius

    This is a simplified virial mass with
    M = sigma_3d**2 * radius / G
    and sigma_3d = sqrt(3) * sigma_1d
    """
    return ((3 * vdisp**2) * radius / constants.G).to(u.M_sun)

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
print("Density from that mass: {0}".format((mvir(0.16*u.pc, np.std(velocities)*u.km/u.s)/(2.8*u.Da)/(4/3.*np.pi*(0.16*u.pc)**3)).to(u.cm**-3)))
velocities_optimistic = ([row['$V_{LSR}$'] for row in tbl22 if row['Object Name'] in cluster_sources_optimistic] +
                         [row['H77a_velocity'] for row in tbl77 if row['ObjectName'] in cluster_sources_optimistic])
print("Mean velocity including weak detections: {0}".format(np.mean(velocities_optimistic)))
print("Velocity dispersion including weak detections: {0}".format(np.std(velocities_optimistic)))
# spanning tree gives 0.09 pc
print("Virial mass in e1e2 using R=0.16 pc (eyeballed): {0}".format(mvir(0.16*u.pc, np.std(velocities_optimistic)*u.km/u.s)))
print("Density from that mass: {0}".format((mvir(0.16*u.pc, np.std(velocities_optimistic)*u.km/u.s)/(2.8*u.Da)/(4/3.*np.pi*(0.16*u.pc)**3)).to(u.cm**-3)))
velocities_optimistic = ([row['$V_{LSR}$'] for row in tbl22 if row['Object Name'] in e1cluster] +
                         [row['H77a_velocity'] for row in tbl77 if row['ObjectName'] in e1cluster])
print("Mean velocity in e1: {0}".format(np.mean(velocities_optimistic)))
print("Velocity dispersion in e1: {0}".format(np.std(velocities_optimistic)))
# spanning tree gives 0.04 pc
print("Virial mass in e1 using R=0.07 pc (eyeballed): {0}".format(mvir(0.07*u.pc, np.std(velocities_optimistic)*u.km/u.s)))
print("Density from that mass: {0}".format((mvir(0.07*u.pc, np.std(velocities_optimistic)*u.km/u.s)/(2.8*u.Da)/(4/3.*np.pi*(0.07*u.pc)**3)).to(u.cm**-3)))

# Are the clusters star or gas dominated?
# see cluster_properties.py
# Implied stellar mass of e1 cluster, assuming 7 stars are OB stars: 695.874032382 +/- 238.642840322 M_sun
# Implied stellar mass of e1 cluster, assuming 7 stars are O stars: 2388.85335601 +/- 752.01450386 M_sun
