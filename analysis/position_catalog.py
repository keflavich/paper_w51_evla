import numpy as np
import paths
import pyregion
from astropy import coordinates
from astropy import units as u
from astropy import table
from astropy.table import Table,Column
from latex_info import latexdict, exp_to_tex, format_float

points = pyregion.open(paths.rpath("pointsource_centroids.reg"))
diffuse = pyregion.open(paths.rpath("diffuse_hii_regions.reg"))

# filter out sources we don't want repeated
diffuse_list = ('d3_diffuse', 'e1', 'e6', 'e7',)
points = [reg for reg in points if reg.attr[1]['text'] not in diffuse_list]

# slower but correct by doing it twice (in case one is fk5 and one is galactic)
coords = coordinates.SkyCoord([coordinates.SkyCoord(*reg.coord_list[:2],
                                                    unit=(u.deg, u.deg),
                                                    frame=reg.coord_format) for
                               reg in points+diffuse])
radii = ([(reg.coord_list[2]) if len(reg.coord_list) > 2 else np.nan
          for reg in points+diffuse]*u.deg).to(u.arcsec)
names = [reg.attr[1]['text'] for reg in points+diffuse]
                              
postbl = Table([Column(data=names, name='Source Name'),
                Column(data=coords.ra.to_string(unit=u.hour, sep=':'), name='RA'),
                Column(data=coords.dec.to_string(unit=u.deg, sep=':'), name='Dec'),
                Column(data=radii, name='Radius'),
                Column(data=(radii*5.1*u.kpc).to(u.pc, u.dimensionless_angles()), name='Phys. Radius'),
               ])
postbl.sort('Source Name')

latexdict['header_start'] = '\label{tab:positions}'
latexdict['caption'] = 'Source Positions'
latexdict['tablefoot'] = '\par\nObjects with name e\#d are the diffuse counterparts to point sources.'
latexdict['col_align'] = 'lllrr'
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
postbl.write(paths.tpath('source_positions.tex'), format='ascii.latex',
             latexdict=latexdict,
             formats={'Radius': lambda x: '-' if np.isnan(x) else '{0:0.2f}'.format(x),
                      'Phys. Radius': lambda x: '-' if np.isnan(x) else '{0:0.2f}'.format(x),
                     },
            )
