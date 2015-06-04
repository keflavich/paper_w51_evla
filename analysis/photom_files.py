from paths import datapath

files = {
        # too low-res '1.4 GHz Epoch 1': 'W51-LBAND-feathered_ABCD.fits',
        '2.5 GHz Epoch 2': datapath+"W51_12B-365_2to3GHz_continuum_uniform.image.fits",
        '3.5 GHz Epoch 2': datapath+"W51_12B-365_3to4GHz_continuum_uniform.image.fits",
        '4.9 GHz Epoch 2': datapath+"W51_12B-365_4.4to5.4GHz_continuum_uniform.image.fits",
        '4.9 GHz Epoch 1': datapath+'W51-CBAND-feathered.fits',
        '4.9 GHz Epoch 3': datapath+'W51Ku_C_Aarray_continuum_2048_low_uniform.clean.image.fits',
        # 4096 A/C seem to suppress large scales WAY too much
        #'4.9 GHz Epoch 3': datapath+'W51C_ACarray_continuum_4096_low_uniform.clean.image.fits',
        '5.9 GHz Epoch 2': datapath+"W51_12B-365_5.4to6.4GHz_continuum_uniform.image.fits",
        '5.9 GHz Epoch 3': datapath+'W51Ku_C_Aarray_continuum_2048_high_uniform.clean.image.fits',
        #'5.9 GHz Epoch 3': datapath+'W51C_ACarray_continuum_4096_high_uniform.clean.image.fits',
        '8.4 GHz Epoch 1': datapath+'W51-X-ABCD-S1.VTESS.VTC.DAVID-MEH.fits',
        '14.1 GHz Epoch 2':datapath+"W51Ku_BDarray_continuum_2048_high_uniform.hires.clean.image.fits",
        #'13.0 GHz':datapath1+"w51_ku/W51Ku_BDarray_continuum_2048_both_uniform_GBTmodel.hires.clean.image.fits",
        '12.6 GHz Epoch 2':datapath+"W51Ku_BDarray_continuum_2048_low_uniform.hires.clean.image.fits",
        '22.5 GHz Epoch 1':datapath+"W51-K-B.S1-ICLN.DAVID-MEH.fits",
        }
