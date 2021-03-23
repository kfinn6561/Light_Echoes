'''
Created on Jun 30, 2015

@author: kieran
'''
'''
Created on May 25, 2015

RUN SOURCE EXPORT.SOURCEME BEFORE STARTING


@author: kieran
'''
from spectral_tools import *
import pylab as p
import numpy as np
from general_tools import progress_bar,pload,pdump
from functions import integrate,integrate_data,log_n

def get_continuum(x,y):
    a,b,c,d=np.polyfit(x,y,3)
    out=a*(x**3)+b*(x**2)+c*x+d
    return out

def scale(x):#scale such that the average value is 1
    x=np.array(x)
    norm=np.sum(x)/len(x)
    return x/norm

def calc_photo(wlength,flux,band):
    f=get_band(band)
    y=flux*f(wlength)
    return integrate_data(wlength, y)

mag_base=np.power(100,1./5)

supernova='sn1993j'
reduction_type='dered-warp'    

print('Loading the data')
sp=spec_plotter()
sp.loadspeclist(supernova)
sp.specreductiontype=reduction_type
sp.loadspectra()


phase,wlength,flux=sp.extract_data()
BmV=[]
VmH=[]
RmI=[]
for i in range(len(phase)):
    progress_bar(i, len(phase))
    photo={}
    for band in ['B','V','H','R','I']:
        photo[band]=calc_photo(wlength[i],flux[i],band)
    BmV.append(-log_n(photo['B']/photo['V'],mag_base))
    VmH.append(-log_n(photo['V']/photo['H'],mag_base))
    RmI.append(-log_n(photo['R']/photo['I'],mag_base))
    
    cont=get_continuum(wlength[i],flux[i])
    out=scale(flux[i]/cont)-1
    pdump([wlength[i],out],'meanspec/meanspectest_spec_%g.sav' %phase[i])

out={'B-V':[phase,BmV,np.zeros_like(phase)],
     'V-H':[phase,VmH,np.zeros_like(phase)],
     'R-I':[phase,RmI,np.zeros_like(phase)]}
pdump(out,'fake_photometry.pkl')
    









