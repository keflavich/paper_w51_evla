"""
Run ptsrc_photom interactively first
and %run -i this
"""
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle
import pylab as pl
import paths
from astropy import log



sorter = lambda x: float(x.split()[0])
freqs = sorted([float(k.split()[0])*u.Unit(k.split()[1]) for k in fluxes.keys()])
#freqs = np.array(sorted(fluxes.keys(),key=sorter)) *u.GHz
sorted_keys = sorted(fluxes.keys(),key=sorter) #[" ".join([str(x.value), x.unit.to_string()]) for x in freqs]
freqs = np.array([f.value for f in freqs]) * u.GHz
errs = np.array([error[k] for k in sorted_keys])

ep1 = np.array(['Epoch 1' in k for k in sorted_keys],dtype='bool')
ep2 = np.array(['Epoch 2' in k for k in sorted_keys],dtype='bool')
ep3 = np.array(['Epoch 3' in k for k in sorted_keys],dtype='bool')

farr = np.linspace(1,40,100)

#for ii,sid in enumerate(string.ascii_uppercase[:7]):
for ii,sid in enumerate(fluxes['2.5 GHz Epoch 2'].keys()):

    log.info("SID: {0}".format(sid))
    fig = pl.figure(ii,figsize=(10,13))
    fig.clf()

    #fplot = [fluxes[k]['A'+sid] for k in sorted_keys]
    fplot = np.array([peaks[k][sid] for k in sorted_keys])
    apfplot = np.array([fluxes[k][sid] for k in sorted_keys])
    fmin = np.array([valleys[k][sid] for k in sorted_keys])
    OK = np.isfinite(fplot)
    nOK = np.count_nonzero(OK)

    if nOK > 12:
        spdim1 = 6
        spdim2 = 4
    elif nOK > 9:
        spdim1 = 5
        spdim2 = 4
    else:
        spdim1 = 5
        spdim2 = 3

    ax = fig.add_subplot(spdim1,1,spdim1-1)
    #fig.suptitle(sid,fontsize=20)

    ep1color = (0.2,1,0.2,0.5) # light green
    ep2color = (0,0.2,1,0.5) # blue
    ep3color = (0.8,0.2,0.2,0.5) # red

    ax.errorbar(np.array(freqs)[ep2 & OK], np.array(fplot)[ep2 & OK]*1e3,
                np.array(errs)[ep2 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep2color, color=ep2color)
    ax.errorbar(np.array(freqs)[ep1 & OK], np.array(fplot)[ep1 & OK]*1e3,
                np.array(errs)[ep1 & OK]*1e3, linestyle='none', marker='s',
                color=ep1color, mec='none', mfc=ep1color)
    ax.errorbar(np.array(freqs)[ep3 & OK], np.array(fplot)[ep3 & OK]*1e3,
                np.array(errs)[ep3 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep3color, color=ep3color)
    ax.plot(np.array(freqs)[ep2 & OK], np.array(fplot-fmin)[ep2 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep2color)
    ax.plot(np.array(freqs)[ep1 & OK], np.array(fplot-fmin)[ep1 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep1color)
    ax.plot(np.array(freqs)[ep3 & OK], np.array(fplot-fmin)[ep3 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep3color)
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ind = np.where(freqs == 12.6*u.GHz)[0][0]
    ax.plot(farr,farr**2/(freqs[ind].value**2)*fplot[ind]*1e3, 'k--')
    ax.plot(farr,farr**1/(freqs[ind].value**1)*fplot[ind]*1e3, 'k:')
    ax.fill_between(freqs.value,-errs*1e3,errs*1e3,color=(1,0.1,0.1,0.5))
    ax.fill_between(freqs.value,-errs*1e3*3,errs*1e3*3,color=(1,0.1,0.1,0.2))
    ax.set_xlim(xlims)
    if np.abs(ylims[0]) > ylims[1]:
        ylims = (-np.abs(ylims[1]), ylims[1])
    ax.set_ylim(ylims)
    #ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Peak Flux Density\n (mJy/beam)")
    ax.set_xticklabels([])

    ax = fig.add_subplot(spdim1,1,spdim1)
    ax.errorbar(np.array(freqs)[ep2 & OK], np.array(fplot)[ep2 & OK]*1e3,
                np.array(errs)[ep2 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep2color, color=ep2color)
    ax.errorbar(np.array(freqs)[ep1 & OK], np.array(fplot)[ep1 & OK]*1e3,
                np.array(errs)[ep1 & OK]*1e3, linestyle='none', marker='s',
                color=ep1color, mec='none', mfc=ep1color)
    ax.errorbar(np.array(freqs)[ep3 & OK], np.array(fplot)[ep3 & OK]*1e3,
                np.array(errs)[ep3 & OK]*1e3, linestyle='none', marker='s',
                markeredgecolor='none', markerfacecolor=ep3color, color=ep3color)
    ax.semilogy(np.array(freqs)[ep2 & OK], np.array(fplot-fmin)[ep2 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep2color)
    ax.semilogy(np.array(freqs)[ep1 & OK], np.array(fplot-fmin)[ep1 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep1color)
    ax.semilogy(np.array(freqs)[ep3 & OK], np.array(fplot-fmin)[ep3 & OK]*1e3,
            linestyle='none', marker='o', markeredgecolor='none',
            markerfacecolor=ep3color)
    ylims = ax.get_ylim()
    ind = np.where(freqs == 12.6*u.GHz)[0][0]
    ax.semilogy(farr,farr**2/(freqs[ind].value**2)*fplot[ind]*1e3, 'k--')
    ax.semilogy(farr,farr**1/(freqs[ind].value**1)*fplot[ind]*1e3, 'k:')
    ax.fill_between(freqs.value,1e-3,errs*1e3,color=(1,0.1,0.1,0.5))
    ax.fill_between(freqs.value,1e-3,errs*1e3*3,color=(1,0.1,0.1,0.2))
    ax.set_xlim(xlims)
    if ylims[0] < 0.05:
        ylims = (0.05, ylims[1])
    ax.set_ylim(ylims)
    ax.set_xlabel("Frequency (GHz)")
    #ax.set_ylabel("Peak Flux Density\n (mJy/beam)")

    ax = fig.add_subplot(spdim1,1,spdim1)


    for jj,(key,isOK) in enumerate(zip(sorted_keys, OK)):

        # skip out-of-field
        if not isOK:
            continue

        sh = cutouts[key][sid].shape

        # skip empties
        if any([s==0 for s in sh]):
            continue

        spnum = np.where(OK)[0].tolist().index(jj)

        ax = fig.add_subplot(spdim1,spdim2,1+spnum)
        rightside = spnum%spdim2 == (spdim2-1) # is this plot on the right edge?
        im = ax.imshow(cutouts[key][sid]*1000, cmap=pl.cm.gray_r) # convert to mJy
        #ax.contour(gfits[k][sid], levels=np.array([2,5,10,50])*errors[k][sid], colors=['b']*10)
        ax.add_artist(Circle((sh[1]/2.,sh[0]/2.),radius=sh[0]/6.,facecolor='none',edgecolor='#FF0000',alpha=0.5))
        ax.set_xticks([])
        ax.set_yticks([])

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cb = pl.colorbar(im, cax=cax)
        if rightside:
            cb.set_label('mJy/beam')
        #ax.annotate(key,[0.9,0.1],xycoords='axes fraction',color='k',ha='right',weight='bold')
        ax.set_title(key.replace("Epoch ","E"))
        #ax.annotate(key.replace("Epoch ","E"), [0.9, 0.1], xycoords='axes fraction', color='k',
        #            ha='right', weight='bold')

    pl.subplots_adjust(hspace=0.1,wspace=0.5)
    fig.savefig(paths.ptsrc_sedpath('%s_SED.pdf' % sid),
                bbox_inches='tight')
    fig.savefig(paths.ptsrc_sedpath('%s_SED.png' % sid),
                bbox_inches='tight')

    #break

pl.show()
