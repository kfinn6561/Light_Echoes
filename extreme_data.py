'''
Created on Jul 22, 2015

@author: kieran
'''
from spectral_tools import spec_plotter,Sav_Template,apply_window,get_band
from functions import log_n
from echo_tools import *
import pylab as p
import os
from general_tools import progress_bar,pload,pdump,reduce_range
import sys

import matplotlib
matplotlib.rcParams.update({'font.size': 16})#increase font size

T=2009.
d0=11000.#distance to Cas a (=3.4kpc
mag_base=100.**(1./5)#base of the magnitude system


def scale(x,*args):#scale such that the average value is 1
    x=np.array(x)
    norm=np.abs(np.sum(x)/len(x))
    out=[x/norm]
    for arg in args:
        out.append(arg/norm)
    return out

class Line():
    def __init__(self,name,lambda_min,lambda_max,l0,band):
        self.name=name
        self.minl=lambda_min
        self.maxl=lambda_max
        self.band=band
        self.l0=l0

#                  sigma,   rho_0,      z0,         alpha,  distance
le_2116=dust_sheet(0.031,   627.9561,   437.11,     9.,     d0,     T,psf=0.97,offset=-0.42,delta=-39.97,profile=gaussian)
le_2521=dust_sheet(0.093,   317.9613,   -9.89,      54.,    d0,     T,psf=1.13,offset=-0.01,delta=-0.91,profile=gaussian)
le_3923=dust_sheet(0.614,   1203.8904,  2045.38,    7.,     d0,     T,psf=1.65,offset=-0.35,delta=-35.36,profile=gaussian)

les=[le_2116,le_2521,le_3923]
colours=['r','b','y','k']#last colour is for unmodified

lines=[Line('H_alpha',5900.,6600.,6564.,'r'),
       Line('He',5300.,6000.,5876.,'V'),
       Line('He',5300.,6000.,5876.,'r'),
       Line('Ca_III',7800.,8800.,8579.,'i')]



bands=['B','V','r','i']
line_colors={'H_alpha':'b','He':'r','O':'g','Ca_III':'c'}
x=np.linspace(3000,10000,1000)

sn=Sav_Template('IIb_1specperSN')
sn.get_lightcurve()
ts_days,wlengths,_=sn.extract_data()#finds the range of phases and wavelthengths required

full_lambdas=[]
for a in wlengths:
    for b in a:
        full_lambdas.append(b)
full_lambdas=np.sort(list(set(full_lambdas)))#remove duplicates and sort

def get_windows(band,recalculate=False):
    fname='%s_band_extreme_window_data.dat' %band
    if recalculate:#forces recalculation
        os.remove(fname)        
    try:
        windows=pload(fname)
        print '%s band window functions read from file' %band
    except IOError:
        print '\n\n'
        print 'calculating window functions for the %s band' %band
        windows=[[] for le in les]
        slit=Slit(0.,1.,subslits=8)
        ts_years=ts_days/365.#todo change echo_tools to work in days
        sn.get_lightcurve(band)
        source=Source(sn.lightcurve[band],min(sn.lc_phases)/365.,max(sn.lc_phases/365.),T-300)
        for i in range(len(les)):
            print '\n'
            print 'Dust sheet %d of %d' %(i+1,len(les))
            for j in range(len(ts_years)):
                progress_bar(j,len(ts_years))
                windows[i].append(les[i].window(ts_years[j],T,slit,source=source))
            windows[i]=normalise(np.array(windows[i]))
        pdump(windows,fname)
        print '\n'
    return windows


print '\n'


print 'producing spectra'

out={}
for line in lines:
    windows=get_windows(line.band)
    test=np.array([len(windows[i])==len(ts_days) for i in range(len(windows))])
    if not test.all():
        print 'window functions out of date, need to recalculate'
        windows=get_windows(line.band,recalculate=True)
    windows.append(np.ones_like(windows[0]))#for unmodified
    lambdas=lambdas=reduce_range(full_lambdas,line.minl,line.maxl)
    
    out[line.name]={'lambdas':lambdas,'l0':line.l0}#for output
    
    spectra=[]
    mins=[]
    maxs=[]
    extreme_min=[]
    extreme_max=[]
    for i in range(len(windows)):
        spectra.append(apply_window(sn,windows[i], ts_days, lambdas,band=line.band))
        mins.append(apply_window(sn,windows[i], ts_days, lambdas,band=line.band,minmax='min'))
        maxs.append(apply_window(sn,windows[i], ts_days, lambdas,band=line.band,minmax='max'))
        extreme_min.append(apply_window(sn,windows[i], ts_days, lambdas,band=line.band,minmax='extreme_min'))
        extreme_max.append(apply_window(sn,windows[i], ts_days, lambdas,band=line.band,minmax='extreme_max'))
    
    labels=['le2116','le2521','le3923','unmodified']
   
    
    for i in range(len(spectra)):
        x,y,z,minn,maxx=scale(spectra[i],mins[i],maxs[i],extreme_min[i],extreme_max[i])
        out[line.name][labels[i]]=[x,y,z,minn,maxx]
pdump(out,'extreme.pkl')
p.show()




