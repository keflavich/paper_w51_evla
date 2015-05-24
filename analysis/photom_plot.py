"""
Run ptsrc_photom interactively first
"""
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle
import pylab as pl
import paths



sorter = lambda x: float(x.split()[0])
freqs = sorted([float(k.split()[0])*u.Unit(k.split()[1]) for k in fluxes.keys()])
#freqs = np.array(sorted(fluxes.keys(),key=sorter)) *u.GHz
sorted_keys = sorted(fluxes.keys(),key=sorter) #[" ".join([str(x.value), x.unit.to_string()]) for x in freqs]
freqs = np.array([f.value for f in freqs]) * u.GHz
errs = np.array([error[k] for k in sorted_keys])

ep1 = np.array(['Epoch 1' in k for k in sorted_keys],dtype='bool')
ep2 = np.array(['Epoch 2' in k for k in sorted_keys],dtype='bool')
ep3 = np.array(['Epoch 3' in k for k in sorted_keys],dtype='bool')

farr = np.linspace(1,25,100)

#for ii,sid in enumerate(string.ascii_uppercase[:7]):
for ii,sid in enumerate(fluxes['2.5 GHz Epoch 2'].keys()):

    fig = pl.figure(ii,figsize=(8.5,12))
    fig.clf()

    #fplot = [fluxes[k]['A'+sid] for k in sorted_keys]
    fplot = np.array([peaks[k][sid] for k in sorted_keys])
    apfplot = np.array([fluxes[k][sid] for k in sorted_keys])
    fmin = np.array([valleys[k][sid] for k in sorted_keys])
    if len(sorted_keys) > 9:
        spdim2 = 4
    else:
        spdim2 = 3

    ax = fig.add_subplot(4,1,4)
    fig.suptitle(sid,fontsize=20)
    ax.errorbar(np.array(freqs)[ep2], np.array(fplot)[ep2]*1e3,
             np.array(errs)[ep2]*1e3, linestyle='none', marker='s')
    ax.errorbar(np.array(freqs)[ep1], np.array(fplot)[ep1]*1e3,
             np.array(errs)[ep1]*1e3, linestyle='none', marker='s',
             color=(0.2,1,0.2,0.5), mec='none', mfc=(0.2,1,0.2,0.5))
    ax.errorbar(np.array(freqs)[ep3], np.array(fplot)[ep3]*1e3,
             np.array(errs)[ep3]*1e3, linestyle='none', marker='s')
    ax.plot(np.array(freqs)[ep2], np.array(fplot-fmin)[ep2]*1e3, linestyle='none',
         marker='o', markeredgecolor='none', markerfacecolor=(0,0.2,1,0.5))
    ax.plot(np.array(freqs)[ep1], np.array(fplot-fmin)[ep1]*1e3, linestyle='none',
         marker='o', markeredgecolor='none', markerfacecolor=(0.1,0.8,0.1,0.5))
    ax.plot(np.array(freqs)[ep3], np.array(fplot-fmin)[ep3]*1e3, linestyle='none',
         marker='o', markeredgecolor='none', markerfacecolor=(0.8,0.2,0.2,0.5))
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ax.plot(farr,farr**2/(freqs[-3].value**2)*fplot[-3]*1e3, 'k--')
    ax.plot(farr,farr**1/(freqs[-3].value**1)*fplot[-3]*1e3, 'k:')
    ax.fill_between(freqs.value,-errs*1e3,errs*1e3,color=(1,0.1,0.1,0.5))
    ax.fill_between(freqs.value,-errs*1e3*3,errs*1e3*3,color=(1,0.1,0.1,0.2))
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xlabel("Frequency (GHz)")
    ax.set_ylabel("Peak Flux Density\n (mJy/beam)")

    for jj,k in enumerate(sorted_keys):
        sh = cutouts[k][sid].shape

        # skip empties
        if any([s==0 for s in sh]):
            continue

        ax = fig.add_subplot(4,spdim2,1+jj)
        im = ax.imshow(cutouts[k][sid]*1000) # convert to mJy
        ax.contour(gfits[k][sid], levels=np.array([2,5,10,50])*errors[k][sid], colors=['w']*10)
        ax.add_artist(Circle((sh[1]/2.,sh[0]/2.),radius=sh[0]/6.,facecolor='none',edgecolor='gray',alpha=0.5))
        ax.set_xticks([])
        ax.set_yticks([])

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        cb = pl.colorbar(im, cax=cax)
        cb.set_label('mJy/beam')
        #ax.annotate(k,[0.9,0.1],xycoords='axes fraction',color='k',ha='right',weight='bold')
        ax.annotate(k,[0.9,0.1],xycoords='axes fraction',color='w',ha='right',weight='bold')

    pl.subplots_adjust(hspace=0.05,wspace=0.15)
    fig.savefig(os.path.join(paths.ptsrc_sedpath,'%s_SED.pdf' % sid),
                bbox_inches='tight')
    fig.savefig(os.path.join(paths.ptsrc_sedpath,'%s_SED.png' % sid),
                bbox_inches='tight')

    #break

pl.show()
