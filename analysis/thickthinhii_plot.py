"""
Run ptsrc_photom interactively first
and %run -i this
...or I'll try the import approach...
"""
from ptsrc_photom import fluxes, error, peaks, valleys
from astropy import units as u
import numpy as np
import pylab as pl
import paths
from astropy import log
import hIIregion


sorter = lambda x: float(x.split()[0])
freqs = sorted([float(k.split()[0])*u.Unit(k.split()[1]) for k in fluxes.keys()])
#freqs = np.array(sorted(fluxes.keys(),key=sorter)) *u.GHz
sorted_keys = sorted(fluxes.keys(),key=sorter) #[" ".join([str(x.value), x.unit.to_string()]) for x in freqs]
freqs = np.array([f.value for f in freqs]) * u.GHz
errs = np.array([error[k] for k in sorted_keys])

ep1 = np.array(['Epoch 1' in k for k in sorted_keys],dtype='bool')
ep2 = np.array(['Epoch 2' in k for k in sorted_keys],dtype='bool')
ep3 = np.array(['Epoch 3' in k for k in sorted_keys],dtype='bool')
ep4 = np.array(['Epoch 4' in k for k in sorted_keys],dtype='bool')

farr = np.linspace(1,250,500)

#for ii,sid in enumerate(string.ascii_uppercase[:7]):
for ii,sid in enumerate(fluxes['2.5 GHz Epoch 2'].keys()):

    log.info("SID: {0}".format(sid))
    fig = pl.figure(1,figsize=(10,13))
    fig.clf()

    #fplot = [fluxes[k]['A'+sid] for k in sorted_keys]
    fplot = np.array([peaks[k][sid] for k in sorted_keys])
    apfplot = np.array([fluxes[k][sid] for k in sorted_keys])
    fmin = np.array([valleys[k][sid] for k in sorted_keys])
    OK = np.isfinite(fplot)
    nOK = np.count_nonzero(OK)


    ax = fig.add_subplot(2,1,1)
    #fig.suptitle(sid,fontsize=20)

    ep1color = (0.2,1,0.2,0.5) # light green
    ep2color = (0,0.2,1,0.5) # blue
    ep3color = (0.8,0.2,0.2,0.5) # red
    ep4color = (0.0,0.0,0.0,0.5) # black

    ax.errorbar(np.array(freqs)[ep2 & OK], np.array(fplot)[ep2 & OK]*1e3,
                np.array(errs)[ep2 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep2color, color=ep2color)
    ax.errorbar(np.array(freqs)[ep1 & OK], np.array(fplot)[ep1 & OK]*1e3,
                np.array(errs)[ep1 & OK]*1e3, linestyle='none', marker='s',
                color=ep1color, mec='none', mfc=ep1color)
    ax.errorbar(np.array(freqs)[ep3 & OK], np.array(fplot)[ep3 & OK]*1e3,
                np.array(errs)[ep3 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep3color, color=ep3color)
    ax.errorbar(np.array(freqs)[ep4 & OK], np.array(fplot)[ep4 & OK]*1e3,
                np.array(errs)[ep4 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep4color, color=ep4color)
    ax.plot(np.array(freqs)[ep2 & OK], np.array(fplot-fmin)[ep2 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep2color)
    ax.plot(np.array(freqs)[ep1 & OK], np.array(fplot-fmin)[ep1 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep1color)
    ax.plot(np.array(freqs)[ep3 & OK], np.array(fplot-fmin)[ep3 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep3color)
    ax.plot(np.array(freqs)[ep4 & OK], np.array(fplot-fmin)[ep4 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep4color)
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ind = np.where(freqs == 12.6*u.GHz)[0][0]
    #ax.plot(farr,farr**2/(freqs[ind].value**2)*fplot[ind]*1e3, 'k--')
    #ax.plot(farr,farr**1/(freqs[ind].value**1)*fplot[ind]*1e3, 'k:')
    ax.plot(farr,farr**2/(226**2)*fplot[-1]*1e3, 'b--')


    mp = hIIregion.mpfit(hIIregion.mpfitfun(np.array(freqs)[OK&(ep2|ep3)],
                                            np.array(fplot)[OK&(ep2|ep3)],
                                            np.array(errs)[OK&(ep2|ep3)]),
                         xall=[7,1.0e-6], quiet=1)
    em = mp.params[0]
    normfac = mp.params[1]
    
    if em < 3:
        em = 4

    electron_temperature = 7500
    hiimodel_bf = hIIregion.inufit(farr, em, normfac)     *1e3
    hiimodel_x10 = hIIregion.inufit(farr, em+1, normfac) *1e3
    hiimodel_d10 = hIIregion.inufit(farr, em-1, normfac) *1e3
    ax.plot(farr, hiimodel_bf, 'k--')
    ax.plot(farr, hiimodel_x10, 'k:')
    ax.plot(farr, hiimodel_d10, 'k-.')

    ax.fill_between(freqs.value,-errs*1e3,errs*1e3,color=(1,0.1,0.1,0.5))
    ax.fill_between(freqs.value,-errs*1e3*3,errs*1e3*3,color=(1,0.1,0.1,0.2))
    ax.set_xlim(xlims)
    if np.abs(ylims[0]) > ylims[1]:
        ylims = (-np.abs(ylims[1]), ylims[1])
    ax.set_ylim(ylims)
    #ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Peak Flux Density\n (mJy/beam)")
    ax.set_xticklabels([])

    ax = fig.add_subplot(2,1,2)
    ax.errorbar(np.array(freqs)[ep2 & OK], np.array(fplot)[ep2 & OK]*1e3,
                np.array(errs)[ep2 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep2color, color=ep2color)
    ax.errorbar(np.array(freqs)[ep1 & OK], np.array(fplot)[ep1 & OK]*1e3,
                np.array(errs)[ep1 & OK]*1e3, linestyle='none', marker='s',
                color=ep1color, mec='none', mfc=ep1color)
    ax.errorbar(np.array(freqs)[ep3 & OK], np.array(fplot)[ep3 & OK]*1e3,
                np.array(errs)[ep3 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep3color, color=ep3color)
    ax.errorbar(np.array(freqs)[ep4 & OK], np.array(fplot)[ep4 & OK]*1e3,
                np.array(errs)[ep4 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep4color, color=ep4color)
    ax.semilogy(np.array(freqs)[ep2 & OK], np.array(fplot-fmin)[ep2 & OK]*1e3,
                linestyle='none', marker='o', markeredgecolor='none',
                markerfacecolor=ep2color)
    ax.semilogy(np.array(freqs)[ep1 & OK], np.array(fplot-fmin)[ep1 & OK]*1e3,
                linestyle='none', marker='o', markeredgecolor='none',
                markerfacecolor=ep1color)
    ax.semilogy(np.array(freqs)[ep3 & OK], np.array(fplot-fmin)[ep3 & OK]*1e3,
                linestyle='none', marker='o', markeredgecolor='none',
                markerfacecolor=ep3color)
    ax.semilogy(np.array(freqs)[ep4 & OK], np.array(fplot-fmin)[ep4 & OK]*1e3,
                linestyle='none', marker='o', markeredgecolor='none',
                markerfacecolor=ep4color)
    ylims = ax.get_ylim()
    ind = np.where(freqs == 12.6*u.GHz)[0][0]
    ax.semilogy(farr,farr**2/(freqs[ind].value**2)*fplot[ind]*1e3, 'k--')
    ax.semilogy(farr,farr**1/(freqs[ind].value**1)*fplot[ind]*1e3, 'k:')
    ax.semilogy(farr,farr**2/(226**2)*fplot[-1]*1e3, 'b--')
    ax.fill_between(freqs.value,1e-3,errs*1e3,color=(1,0.1,0.1,0.5))
    ax.fill_between(freqs.value,1e-3,errs*1e3*3,color=(1,0.1,0.1,0.2))
    ax.set_xlim(xlims)
    if ylims[0] < 0.05:
        ylims = (0.05, ylims[1])
    ax.set_ylim(ylims)
    ax.set_xlabel("Frequency (GHz)")
    #ax.set_ylabel("Peak Flux Density\n (mJy/beam)")

    fig.savefig(paths.ptsrc_sedpath('%s_SED_HII.png' % sid),
                bbox_inches='tight')
