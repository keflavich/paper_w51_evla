import numpy as np
import paths
from astropy import table
from astropy import units as u


tbl22 = table.Table.read(paths.tpath('H2CO22_emission_spectral_fits.ecsv'), format='ascii.ecsv')
tbl77 = table.Table.read(paths.tpath('H77a_spectral_fits.ipac'), format='ascii.ipac')

cluster_sources = ['e8mol', 'e10mol', 'e1', 'e2', 'e3', 'e4', 'e8']
cluster_sources_optimistic = ['e8mol', 'e10mol', 'e1', 'e2', 'e3', 'e4', 'e8', 'e9', 'e10']
e1cluster = ['e8mol', 'e10mol', 'e1', 'e3', 'e4', 'e8',]

velocities = ([row['$V_{LSR}$'] for row in tbl22 if row['Object Name'] in cluster_sources] +
              [row['H77a_velocity'] for row in tbl77 if row['ObjectName'] in cluster_sources])
print("Mean velocity including only solid detections: {0}".format(np.mean(velocities)))
print("Velocity dispersion including only solid detections: {0}".format(np.std(velocities)))
velocities_optimistic = ([row['$V_{LSR}$'] for row in tbl22 if row['Object Name'] in cluster_sources_optimistic] +
                         [row['H77a_velocity'] for row in tbl77 if row['ObjectName'] in cluster_sources_optimistic])
print("Mean velocity including weak detections: {0}".format(np.mean(velocities_optimistic)))
print("Velocity dispersion including weak detections: {0}".format(np.std(velocities_optimistic)))
velocities_optimistic = ([row['$V_{LSR}$'] for row in tbl22 if row['Object Name'] in e1cluster] +
                         [row['H77a_velocity'] for row in tbl77 if row['ObjectName'] in e1cluster])
print("Mean velocity in e1: {0}".format(np.mean(velocities_optimistic)))
print("Velocity dispersion in e1: {0}".format(np.std(velocities_optimistic)))
