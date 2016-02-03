from __future__ import print_function
import pyregion
import photutils
import astropy.wcs
from astropy.io import fits
import numpy as np
from FITS_tools.strip_headers import flatten_header
from astropy import units as u
import string
from astropy import coordinates
from astropy import table
import gaussfitter
import paths
from paths import datapath
import radio_beam
from astropy import time
from astropy.utils.console import ProgressBar
from rounded import rounded
from latex_info import latexdict, exp_to_tex, format_float
from photom_files import files

reglist = pyregion.open(paths.rpath('pointsource_centroids.reg'))

max_rad = 4 # arcsec

fluxes = {}
peaks = {}
valleys = {}
cutouts = {}
error = {}
mpfits = {}
gpars = {}
gparerrs = {}
gfits = {}
obsdate = {}
beams = {}
frequencies = {}
reg_centers = {}

for freq,fn in files.items():
    fluxes[freq] = {}
    peaks[freq] = {}
    valleys[freq] = {}
    cutouts[freq] = {}
    mpfits[freq] = {}
    gpars[freq] = {}
    gparerrs[freq] = {}
    gfits[freq] = {}
    reg_centers[freq] = {}
    
    data = fits.getdata(fn).squeeze()
    header = flatten_header(fits.getheader(fn))
    wcs = astropy.wcs.WCS(header).sub([astropy.wcs.WCSSUB_CELESTIAL])
    beam = radio_beam.Beam.from_fits_header(header)
    beams[freq] = beam
    frequencies[freq] = header.get('CRVAL3') or header.get('ACRVAL3')

    if wcs.wcs.radesys == 'FK4' or wcs.wcs.equinox == 1950:
        (blx,bly),(tr_x,tr_y) = wcs.wcs_world2pix([(290.35090,14.42627),(290.34560,14.43126),],0)
    else:
        (blx,bly),(tr_x,tr_y) = wcs.wcs_world2pix([(290.92445,14.524376),(290.91912,14.529338)],0)
    error[freq] = data[bly:tr_y,blx:tr_x].std()
    obsdate[freq] = fits.getheader(fn)['DATE-OBS']
    print("file: {0}".format(fn))

    for reg in ProgressBar(reglist):
        if reg.name == 'circle':
            ra,dec,rad = reg.coord_list
        elif reg.name == 'point':
            ra,dec = reg.coord_list
            rad = beam.major.to(u.deg).value

        if wcs.wcs.radesys == 'FK4' or wcs.wcs.equinox == 1950:
            C = coordinates.SkyCoord(ra,dec,unit=(u.deg,u.deg), frame='icrs').fk4
            ra,dec = C.fk4.ra.deg,C.fk4.dec.deg

        rd = [ra,dec] + [0]*(wcs.wcs.naxis-2)
        xc,yc = wcs.wcs_world2pix([rd],0)[0][:2]
        #pixscale = np.abs(wcs.wcs.get_cdelt()[0])
        pixscale = (wcs.pixel_scale_matrix.diagonal()**2).sum()**0.5
        rp = rad / pixscale  # radius in pixels
        aperture = photutils.CircularAperture(positions=[xc, yc], r=rp)
        aperture_data = photutils.aperture_photometry(data=data, apertures=aperture)
        flux = aperture_data[0]['aperture_sum']
        #flux = np.array(photutils.aperture_photometry(data=data,
        #                                              positions=[xc, yc],
        #                                              apertures=('circular',
        #                                                         rp))[0]['aperture_sum'])

        name = reg.attr[1]['text']
        fluxes[freq][name] = flux / (np.pi * rp**2)
        # co: short for cutout
        co = data[int(yc-3*rp):int(yc+3*rp+1),
                  int(xc-3*rp):int(xc+3*rp+1)]
        max_rp = max_rad/3600./pixscale
        cotoplot = data[int(yc-max_rp):int(yc+max_rp+1),
                        int(xc-max_rp):int(xc+max_rp+1)]
        cutouts[freq][name] = (co, cotoplot, pixscale)
        shiftx = xc-int(xc-max_rp)-(xc-int(xc-3*rp))
        shifty = yc-int(yc-max_rp)-(yc-int(yc-3*rp))
        reg_centers[freq][name] = (ra, dec, xc-int(xc-3*rp), yc-int(yc-3*rp),
                                   rad, rp, shiftx, shifty)

        yy,xx = np.indices(co.shape)
        # rr = radius grid
        rr = ((xx-rp*3)**2+(yy-rp*3)**2)**0.5
        # mask = pixels within 1 beam radius
        mask = rr<rp
        if mask.sum():
            peaks[freq][name] = co[mask].max()
            valleys[freq][name] = co.min()
            params = [co.min(), co[mask].max(), 3*rp, 3*rp, 2.0, 2.0, 45.0]
            mp, fitimg = gaussfitter.gaussfit(co, params=params,
                                              err=error[freq],
                                              returnmp=True, rotate=True,
                                              vheight=True, circle=False,
                                              returnfitimage=True)
            mpfits[freq][name] = mp
            gpars[freq][name] = mp.params
            gparerrs[freq][name] = mp.perror
            gfits[freq][name] = fitimg

            wcsX,wcsY = wcs.wcs_pix2world(mp.params[2]+int(xc-3*rp), mp.params[3]+int(yc-3*rp), 0)
            center_coord = coordinates.SkyCoord(wcsX*u.deg, wcsY*u.deg,
                                                frame=wcs.wcs.radesys.lower())
            gpars[freq][name][2] = center_coord.transform_to('fk5').ra.deg
            gpars[freq][name][3] = center_coord.transform_to('fk5').dec.deg
            gpars[freq][name][4] *= pixscale
            gpars[freq][name][5] *= pixscale

            gparerrs[freq][name][2] *= pixscale
            gparerrs[freq][name][3] *= pixscale
            gparerrs[freq][name][4] *= pixscale
            gparerrs[freq][name][5] *= pixscale

        else:
            peaks[freq][name] = np.nan
            valleys[freq][name] = np.nan

            mpfits[freq][name] = ()
            gpars[freq][name] = np.array([np.nan]*7)
            gparerrs[freq][name] = np.array([np.nan]*7)


# hack to make errors like peaks
errors = {freq:
          {name: error[freq] for name in fluxes[freq]}
          for freq in fluxes}
obsdates = {freq:
            {name: obsdate[freq] for name in fluxes[freq]}
            for freq in fluxes}

colname_mappings = {'aperture_flux': fluxes, 'peak_flux': peaks,
                    'cutout_min_flux': valleys, 'local_rms_noise': errors, }
 
tbl = table.Table()
col = table.Column(data=[name for freq in fluxes for name in fluxes[freq]],
                   name='SourceName')
tbl.add_column(col)
col = table.Column(data=[freq for freq in fluxes for name in fluxes[freq]],
                   name='FrequencyName')
tbl.add_column(col)
col = table.Column(data=[obsdates[freq][name] for freq in obsdates for name in obsdates[freq]],
                   name='ObservationDate')
tbl.add_column(col)
col = table.Column(data=[float(freq.split()[0]) for freq in fluxes for name in fluxes[freq]],
                   name='Frequency', unit=u.GHz)
tbl.add_column(col)
col = table.Column(data=[freq.split()[-1] for freq in fluxes for name in fluxes[freq]],
                   name='Epoch')
tbl.add_column(col)
col = table.Column(data=[beams[freq].major.to(u.arcsec).value for freq in beams for name in fluxes[freq]],
                   name='BMAJ', unit=u.arcsec)
tbl.add_column(col)
col = table.Column(data=[beams[freq].minor.to(u.arcsec).value for freq in beams for name in fluxes[freq]],
                   name='BMIN', unit=u.arcsec)
tbl.add_column(col)

for column_name in colname_mappings:
    datadict = colname_mappings[column_name]
    data = [datadict[freq][name]
            for freq in datadict
            for name in datadict[freq]]
    col = table.Column(data=data, unit=u.Jy/u.beam,
                       name=column_name)
    tbl.add_column(col)

gpardict = [('background',u.Jy/u.beam),
            ('amplitude',u.Jy/u.beam),
            ('racen',u.deg),
            ('deccen',u.deg),
            ('xwidth',u.deg),
            ('ywidth',u.deg),
            ('positionangle',u.deg)]
for ii,(parname,unit) in enumerate(gpardict):
    data = [gpars[freq][name][ii]
            for freq in gpars
            for name in gpars[freq]]
    col = table.Column(name='g'+parname, data=data, unit=unit)
    tbl.add_column(col)

data = [mpfits[freq][name].chi2 
        if hasattr(mpfits[freq][name],'chi2') else np.nan
        for freq in mpfits
        for name in mpfits[freq]]
col = table.Column(name='gfit_chi2', data=data)
tbl.add_column(col)

data = [mpfits[freq][name].chi2n
        if hasattr(mpfits[freq][name],'chi2') else np.nan
        for freq in mpfits
        for name in mpfits[freq]]
col = table.Column(name='gfit_chi2_reduced', data=data)
tbl.add_column(col)

rawtbl = tbl.copy()

bad = ((rawtbl['peak_flux']-rawtbl['cutout_min_flux'] < 0) # negative signal = bad
       | ((rawtbl['peak_flux']-rawtbl['cutout_min_flux'])/rawtbl['local_rms_noise'] < 4)
       | np.isnan(rawtbl['peak_flux'])
       | np.isnan(rawtbl['cutout_min_flux'])
      )

# uplims: 'peak_flux', 'aperture_flux', 'local_rms_noise', 'cutout_min_flux',
for colname in ['gbackground', 'gamplitude', 'gracen', 'gdeccen', 'gxwidth',
                'gywidth', 'gpositionangle', 'gfit_chi2', 'gfit_chi2_reduced']:
    tbl[bad][colname] = np.nan

tbl.sort(['SourceName', 'Frequency'])
tbl.write(paths.tpath('EVLA_VLA_PointSourcePhotometry.ipac'), format='ascii.ipac')

## Create latex table

tbl.add_column(table.Column(data=tbl['peak_flux']-tbl['cutout_min_flux'],
                            name='peak_m_background',
                            unit=tbl['peak_flux'].unit))
nondetections = tbl['peak_m_background'] < tbl['local_rms_noise']*3

cols_order = ['SourceName', 'Epoch', 'ObservationDate', 'peak_flux', 'peak_m_background', 'local_rms_noise', 'Frequency']
cols = {'SourceName': 'Object',
        #'FrequencyName': 'Band',
        'ObservationDate': 'Obs. Date',
        'peak_flux': 'Peak $S_{\\nu}$',
        'peak_m_background': 'Peak - Background',
        'local_rms_noise': 'RMS',
       }

textbl = tbl.copy()[cols_order]
textbl.sort(['SourceName', 'Frequency'])
textbl[nondetections]['peak_flux'] = np.nan
textbl[nondetections]['peak_m_background'] = np.nan
textbl['peak_flux'] = ((list(map(lambda x,y: rounded(x,y,extra=0)[0],
                                 textbl['peak_flux'].to(u.mJy/u.beam).value,
                                 textbl['local_rms_noise'].to(u.mJy/u.beam).value))))
textbl['peak_m_background'] = ((list(map(lambda x,y: rounded(x,y,extra=0)[0],
                                         textbl['peak_m_background'].to(u.mJy/u.beam).value,
                                         textbl['local_rms_noise'].to(u.mJy/u.beam).value))))
textbl['local_rms_noise'] = ((list(map(lambda x,y: rounded(x,y,extra=0)[0],
                                       textbl['local_rms_noise'].to(u.mJy/u.beam).value,
                                       textbl['local_rms_noise'].to(u.mJy/u.beam).value))))
for name in ('peak_flux', 'peak_m_background', 'local_rms_noise'):
    textbl[name].unit = u.mJy/u.beam
textbl['SourceName'] = list(map(lambda x: x.replace("_","-"), textbl['SourceName']))

for old,new in cols.items():
    textbl.rename_column(old, new)

latexdict['header_start'] = '\label{tab:contsrcs}'
latexdict['caption'] = 'Continuum Point Sources (excerpt)'
latexdict['tablefoot'] = '\par\nAn excerpt from the point source catalog.  For the full catalog, see Table \\ref{tbl:contsrcs_full}'
textbl[::10].write(paths.tpath('pointsource_photometry.tex'), format='ascii.latex', latexdict=latexdict,
             formats={'RMS': format_float,
                      'Peak $S_{\\nu}$': format_float,
                      'Peak - Background': format_float,
                      'Obs. Date': lambda x: time.Time(x).iso[:10],
                     })

latexdict['header_start'] = '\label{tab:contsrcs_full}'
latexdict['caption'] = 'Continuum Point Sources'
latexdict['tablefoot'] = ''
#latexdict['tabletype'] = 'longtable'
latexdict['tabulartype'] = 'longtable'
textbl.write(paths.tpath('pointsource_photometry_full.tex'), format='ascii.latex', latexdict=latexdict,
             formats={'RMS': format_float,
                      'Peak $S_{\\nu}$': format_float,
                      'Peak - Background': format_float,
                      'Obs. Date': lambda x: time.Time(x).iso[:10],
                     })
