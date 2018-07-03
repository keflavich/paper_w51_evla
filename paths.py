import os
import socket

if socket.gethostname() == 'cleese':
    root = '/scratch/aginsbur/w51/'
else:
    root = os.path.expanduser('~/work/w51')

datapath = os.path.join(root, 'paper_w51_evla/data/')
texpath = os.path.join(root, 'paper_w51_evla/tex/')
figurepath = os.path.join(root, 'paper_w51_evla/figures/')
regpath = os.path.join(root, 'paper_w51_evla/regions/')
gridpath = os.path.join(root, '../h2co/radex/...todo.../')
analysispath = os.path.join(root, 'paper_w51_evla/analysis/')
plotcodepath = os.path.join(root, 'paper_w51_evla/plot_codes/')
observingpath = os.path.join(root, 'paper_w51_evla/observing/')
tablepath = os.path.join(root, 'paper_w51_evla/tables/')
ptsrc_sed_figure_path = os.path.join(figurepath, 'pointsource_seds')
almapath = os.path.join(root, 'alma')

def gpath(x, gridpath=gridpath):
    return os.path.join(gridpath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)

def rpath(x, regpath=regpath):
    return os.path.join(regpath, x)

def opath(x, observingpath=observingpath):
    return os.path.join(observingpath, x)

def pcpath(x, plotcodepath=plotcodepath):
    return os.path.join(plotcodepath, x)

def apath(x, analysispath=analysispath):
    return os.path.join(analysispath, x)

def dpath(x, datapath=datapath):
    return os.path.join(datapath, x)

def dppath(x, datapath=datapath):
    return os.path.join(datapath, 'projections', x)

def tpath(x, tablepath=tablepath):
    return os.path.join(tablepath, x)

def ptsrc_sedpath(x, ptsrc_sed_figure_path=ptsrc_sed_figure_path):
    return os.path.join(ptsrc_sed_figure_path, x)

def texpath(x, texpath=texpath):
    return os.path.join(texpath, x)

def almadpath(x, almapath=almapath):
    return os.path.join(almapath, 'FITS', x)
