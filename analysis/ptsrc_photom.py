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

reglist = pyregion.open('/Users/adam/work/w51/pointish_sources.reg')

datapath1 = '/Volumes/128gbdisk/'
datapath2 = '/Users/adam/work/w51/'

files = {
        # too low-res '1.4 GHz Epoch 1': 'W51-LBAND-feathered_ABCD.fits',
        '2.5 GHz Epoch 2': datapath1+"12B-365/W51_12B-365_2to3GHz_continuum_uniform.image.fits",
        '3.5 GHz Epoch 2': datapath1+"12B-365/W51_12B-365_3to4GHz_continuum_uniform.image.fits",
        '4.9 GHz Epoch 2': datapath1+"12B-365/W51_12B-365_4.4to5.4GHz_continuum_uniform.image.fits",
        '4.9 GHz Epoch 1': datapath2+'W51-CBAND-feathered.fits',
        '4.9 GHz Epoch 3': datapath1+'w51_c_a/W51Ku_C_Aarray_continuum_2048_low_uniform.clean.image.fits',
        '5.9 GHz Epoch 2': datapath1+"12B-365/W51_12B-365_5.4to6.4GHz_continuum_uniform.image.fits",
        '5.9 GHz Epoch 3': datapath1+'w51_c_a/W51Ku_C_Aarray_continuum_2048_high_uniform.clean.image.fits',
        '8.4 GHz Epoch 1': datapath2+'W51-X-ABCD-S1.VTESS.VTC.DAVID-MEH.fits',
        '14.1 GHz Epoch 2':datapath1+"w51_ku/W51Ku_BDarray_continuum_2048_high_uniform.hires.clean.image.fits",
        #'13.0 GHz':datapath1+"w51_ku/W51Ku_BDarray_continuum_2048_both_uniform_GBTmodel.hires.clean.image.fits",
        '12.6 GHz Epoch 2':datapath1+"w51_ku/W51Ku_BDarray_continuum_2048_low_uniform.hires.clean.image.fits",
        '22.5 GHz Epoch 1':datapath2+"W51-K-B.S1-ICLN.DAVID-MEH.fits",
        }


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

for freq,fn in files.iteritems():
    fluxes[freq] = {}
    peaks[freq] = {}
    valleys[freq] = {}
    cutouts[freq] = {}
    mpfits[freq] = {}
    gpars[freq] = {}
    gparerrs[freq] = {}
    gfits[freq] = {}
    
    data = fits.getdata(fn).squeeze()
    wcs = astropy.wcs.WCS(flatten_header(fits.getheader(fn)))

    if wcs.wcs.radesys == 'FK4' or wcs.wcs.equinox == 1950:
        (blx,bly),(tr_x,tr_y) = wcs.wcs_world2pix([(290.35090,14.42627),(290.34560,14.43126),],0)
    else:
        (blx,bly),(tr_x,tr_y) = wcs.wcs_world2pix([(290.92445,14.524376),(290.91912,14.529338)],0)
    error[freq] = data[bly:tr_y,blx:tr_x].std()
    obsdate[freq] = fits.getheader(fn)['DATE-OBS']

    for reg in reglist:
        if reg.name == 'circle':
            ra,dec,rad = reg.coord_list
            if wcs.wcs.radesys == 'FK4' or wcs.wcs.equinox == 1950:
                C = coordinates.ICRS(ra,dec,unit=(u.deg,u.deg)).fk4
                ra,dec = C.fk4.ra.deg,C.fk4.dec.deg
            rd = [ra,dec] + [0]*(wcs.wcs.naxis-2)
            xc,yc = wcs.wcs_world2pix([rd],0)[0][:2]
            pixscale = np.abs(wcs.wcs.get_cdelt()[0])
            rp = rad / pixscale
            #aperture = photutils.CircularAperture(positions=[xc, yc], r=rp)
            #flux = photutils.aperture_photometry(data=data, apertures=aperture)
            flux = np.array(photutils.aperture_photometry(data=data,
                                                          positions=[xc, yc],
                                                          apertures=('circular',
                                                                     rp))[0]['aperture_sum'])

            name = reg.attr[1]['text']
            fluxes[freq][name] = flux / (np.pi * rp**2)
            co = data[int(yc-3*rp):int(yc+3*rp+1),
                      int(xc-3*rp):int(xc+3*rp+1)]
            cutouts[freq][name] = co

            yy,xx = np.indices(co.shape)
            rr = ((xx-rp*3)**2+(yy-rp*3)**2)**0.5
            #print co.shape,rr.shape
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

for column_name in colname_mappings:
    datadict = colname_mappings[column_name]
    data = [datadict[freq][name]
            for freq in datadict
            for name in datadict[freq]]
    col = table.Column(data=data, unit=u.mJy/u.beam,
                       name=column_name)
    tbl.add_column(col)

gpardict = [('background',u.mJy/u.beam),
            ('amplitude',u.mJy/u.beam),
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

tbl.write('EVLA_VLA_PointSourcePhotometry.ipac', format='ascii.ipac')
