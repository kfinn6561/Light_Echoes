'''
Created on Aug 25, 2015

@author: kieran
'''
from spectral_tools import smooth_spectrum
import numpy as np
import pylab as p
from matplotlib import gridspec
from general_tools import progress_bar


def plot_difference(lambdas,s1,s2):
    fig=p.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    p.setp(ax0.get_xticklabels(), visible=False)
    ax0.plot(lambdas,s1)
    ax0.plot(lambdas,s2)
    ax1.plot(lambdas,(s1-s2),'rx')

raw=np.loadtxt('../reducedspec/tyc2116-20090922.006-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
var=np.loadtxt('../reducedspec/tyc2116-20090922.006-variance-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
smooth=np.loadtxt('../reducedspec/tyc2116-20090922.006-smooth-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})

def get_error(vexp):
    tsmooth=smooth_spectrum(raw['l'], raw['f'], var['f'], vexp=vexp)
    out=(tsmooth-smooth['f'])**2
    return np.sum(out)

def make_plot(vexp):
    tsmooth=smooth_spectrum(raw['l'], raw['f'], var['f'], vexp=vexp)
    plot_difference(raw['l'],smooth['f'],tsmooth)


x=np.linspace(0.00333,0.00334,100)
out=[]

for i in range(len(x)):
    progress_bar(i, len(x))
    out.append(get_error(x[i]))
    
p.plot(x,out)
p.show()