import os
import socket

if socket.gethostname() == 'cleese':
    root = '/scratch/aginsbur/apex/'
else:
    root = os.path.expanduser('~/work/w51')

datapath = os.path.join(root, 'paper_w51_evla/data/')
figurepath = os.path.join(root, 'paper_w51_evla/figures/')
regpath = os.path.join(root, 'paper_w51_evla/regions/')
gridpath = os.path.join(root, '../h2co/radex/...todo.../')
analysispath = os.path.join(root, 'paper_w51_evla/analysis/')
plotcodepath = os.path.join(root, 'paper_w51_evla/plot_codes/')
observingpath = os.path.join(root, 'paper_w51_evla/observing/')

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
