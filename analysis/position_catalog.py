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
                                                    frame=reg.coord_format).fk5 for
                               reg in points+diffuse])
radii = ([(reg.coord_list[2]) if len(reg.coord_list) > 2 else np.nan
          for reg in points+diffuse]*u.deg).to(u.arcsec)
names = [reg.attr[1]['text'] for reg in points+diffuse]

with open(paths.tpath("SED_class")) as f:
    lines = f.readlines()
    end = lines.index('\n')
    split2 = lines[0].find('Classification')
    split1 = lines[0].find('SED Class')
    seddata = {x[:split1].strip():"${0}$".format(x[split1:split2].strip()) for x in lines[1:end]}
    SEDclasscolumn = Column(data=[seddata[n] if n in seddata else '$E$' for n in names], name="SED Class")
    classification = {x[:split1].strip():"{0}".format(x[split2:].strip()) for x in lines[1:end]}
    classificationcolumn = Column(data=[classification[n] if n in classification else '-' for n in names], name="Classification")
    footnotes_start = lines.index('Classification Key\n') + 1
    footer = "".join(["{0}{1} \\\\\n"
                      .format("${0}$".format(x[0])
                              if len(x) > 1 and x[1]==':'
                              else "\\newline" if len(x) <= 1
                              else x[0],
                              x[1:].strip()) for x in
                      lines[footnotes_start:]])

postbl = Table([Column(data=names, name='Source Name'),
                Column(data=coords.ra.to_string(unit=u.hour, sep=':', precision=2), name='RA'),
                Column(data=coords.dec.to_string(unit=u.deg, sep=':', precision=1), name='Dec'),
                Column(data=radii, name='Radius'),
                Column(data=(radii*5.1*u.kpc).to(u.pc, u.dimensionless_angles()), name='Phys. Radius'),
                SEDclasscolumn,
                classificationcolumn,
               ])

import natsort
postbl = postbl[natsort.index_natsorted(postbl['Source Name'])]

latexdict['header_start'] = '\label{tab:positions}'
latexdict['caption'] = 'Source Positions'
latexdict['tablefoot'] = ('\par\nObjects with name e\#d are the diffuse '
                          'counterparts to point sources.  '
                          'The absolute positional accuracy is '
                          '$\sim0.2\\arcsec$.  '
                          'Sources with no radius are unresolved, '
                          'with '
                          'upper limits of 0.3\\arcsec (0.007 pc).\\\\\n' +
                         footer)
latexdict['col_align'] = 'lllrrll'
#latexdict['tabletype'] = 'longtable'
#latexdict['tabulartype'] = 'longtable'
postbl.write(paths.tpath('source_positions.tex'), format='ascii.latex',
             latexdict=latexdict,
             formats={'Radius': lambda x: '-' if np.isnan(x) else '{0:0.1f}'.format(x),
                      'Phys. Radius': lambda x: '-' if np.isnan(x) else '{0:0.2f}'.format(x),
                     },
            )
