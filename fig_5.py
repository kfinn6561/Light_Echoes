import pylab as p
import numpy as np
from general_tools import pload,reduce_range,pdump
from scipy.interpolate import interp1d
from copy import copy

import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})#increase font size
mpl.rcParams['figure.max_open_warning']=False
mpl.rcParams['xtick.major.size'] = 20
mpl.rcParams['xtick.major.width'] = 4
mpl.rcParams['xtick.minor.size'] = 10
mpl.rcParams['xtick.minor.width'] = 2#tick size
mpl.rcParams['ytick.major.size'] = 20
mpl.rcParams['ytick.major.width'] = 4
mpl.rcParams['ytick.minor.size'] = 10
mpl.rcParams['ytick.minor.width'] = 2#tick size


Ns=[0.001,0.1,1.,10.,100.]
ls=[0.2,2.,20.,100.,500.]


def scale(x,*args):#scale such that the average maximum is 1, er1 and er2 are errors that are scaled the same way
    x=np.array(x)
    offset=0.
    norm=np.max(x)
    out=[(x-offset)/norm]
    for arg in args:
        out.append((arg-offset)/norm)
    return out

def bin(x,y,N=10):
    outx=[]
    outy=[]
    i=1
    while i*N<=len(x):
        outx.append(np.mean(x[(i-1)*N:i*N]))
        outy.append(np.mean(y[(i-1)*N:i*N]))
        i+=1
    start=(i-1)*N
    if start<len(x):
        outx.append(np.mean(x[start:]))
        outy.append(np.mean(y[start:]))
    return (np.array(outx),np.array(outy))

casa2116=np.loadtxt("../spectra/casales/casa2116-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
casa2521=np.loadtxt("../spectra/casales/casa2521-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
casa3923=np.loadtxt("../spectra/casales/casa3923-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})

raw_2116=np.loadtxt('../reducedspec/tyc2116-20090922.006-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
raw_2521=np.loadtxt('../reducedspec/tyc2521-20090922.005-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
raw_3923=np.loadtxt('../reducedspec/tyc3923-20091023.503-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})

raw_data={'le2116':raw_2116,'le2521':raw_2521,'le3923':raw_3923}

casa_data={'le2116':casa2116,'le2521':casa2521,'le3923':casa3923}
les=['le2116','le2521','le3923']


le_colours={'le2116': 'r','le2521':'b','le3923':'orange'} 

offset=1.
for i in range(len(les)):
    le=les[i]
    rx=raw_data[le]['l']
    ry=raw_data[le]['f']
    
    sx=casa_data[le]['l']
    sy=casa_data[le]['f']
    sv=casa_data[le]['var']
    
    
    sy,ry,sv=scale(sy,ry,sv)
    rx,ry=bin(rx,ry)
    
    ry+=i*offset
    sy+=i*offset
    
    p.plot(rx,ry,'k')
    p.plot(sx,sy,le_colours[le])
    #p.fill_between(sx,sy-sv,sy+sv,color=le_colours[le],alpha=0.5)

p.xlabel(r'Rest Wavelength ($\AA$)')
p.ylabel(r'Scaled $f_\lambda$ + Constant')
p.ylim(-1,4)
p.show()
    
    
    