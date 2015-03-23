import os
import socket

if socket.gethostname() == 'cleese':
    root = '/scratch/aginsbur/w51/'
else:
    root = os.path.expanduser('~/work/w51')

datapath = os.path.join(root, 'paper_w51_evla/data/')
figurepath = os.path.join(root, 'paper_w51_evla/figures/')
regpath = os.path.join(root, 'paper_w51_evla/regions/')
gridpath = os.path.join(root, '../h2co/radex/...todo.../')
analysispath = os.path.join(root, 'paper_w51_evla/analysis/')
plotcodepath = os.path.join(root, 'paper_w51_evla/plot_codes/')
observingpath = os.path.join(root, 'paper_w51_evla/observing/')
tablepath = os.path.join(root, 'paper_w51_evla/tables/')

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
