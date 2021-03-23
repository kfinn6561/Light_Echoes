'''
Created on Sep 16, 2015

@author: kieran
'''
from spectral_tools import *
from functions import log_n
from echo_tools_rect import *
import pylab as p
import numpy as np
import os
from general_tools import progress_bar,pload,pdump,reduce_range
import sys

import matplotlib
matplotlib.rcParams.update({'font.size': 16})#increase font size

T=2009.*365.
d0=11000.#distance to Cas a (=3.4kpc
mag_base=100.**(1./5)#base of the magnitude system
lambda_bin=9.#bins the output by this amount, by skipping values
lambda_bin=3.

if 'data' in sys.argv:
    PLOTTING=False
else:
    PLOTTING=True
    
colours=['r','b','y','k']#last colour is for unmodified

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

def skip_value(lst,x):
    out=[]
    i=0
    while i<len(lst):
        out.append(True)
        x0=lst[i]
        i+=1
        while i<len(lst) and lst[i]-x0<x:
            out.append(False)
            i+=1
    return np.array(out,dtype=bool)
    
        

#                  sigma,   rho_0,      z0,         alpha,  distance
le_2116=dust_sheet(0.031,   627.9561,   437.11,     9.,     d0,     psf=0.97,delta=-39.97,offset=-0.42,profile=gaussian)#delta=-39.97,
le_2521=dust_sheet(0.093,   317.9613,   -9.89,      54.,    d0,     psf=1.13,delta=-0.91,offset=-0.01,profile=gaussian)#delta=-0.91,
le_3923=dust_sheet(0.614,   1203.8904,  2045.38,    7.,     d0,     psf=1.65,delta=-35.36,offset=-0.35,profile=gaussian)#delta=-35.36,

les=[le_2116,le_2521,le_3923]

'''
H_alpha=Line('H_alpha',5900.,6600.,6564.,'r')
He=Line('He',5300.,6000.,5876.,'V')
Ca=Line('Ca_III',7800.,8800.,8579.,'i')
lines=[H_alpha,He,Ca]
line_plotnames={'He':'He 5876','H_alpha':r'$H\alpha$','Ca_III':'Ca II'}'''
lines=[Line('H_alpha',5900.,6600.,6564.,'r'),
       Line('He',5300.,6000.,5876.,'V'),
       Line('Ca_III',7800.,8800.,8579.,'i')]
#Line('Fe_II',4900.,5300.,5169,'V')

line_plotnames={'He':'He 5876','H_alpha':r'$H\alpha$','Ca_III':'Ca II','Fe_II':'Fe 5169'}

band_dictionary={'r':'R_BVRI','i':'I_BVRI'}

sn2008ax=flm_Template('sn2008ax',extrapolate_lightcurve=500.)
sn1996cb=flm_Template('sn1996cb')
sn1993J=flm_Template('sn1993J',extrapolate_lightcurve=500.)
sn2011dh=flm_Template('sn2011dh',extrapolate_lightcurve=500.)
sn2011ei=flm_Template('sn2011ei',extrapolate_lightcurve=500.)
sn2003bg=flm_Template('sn2003bg')

#sn2011ei_without=flm_Template('sn2011ei_no33',extrapolate_lightcurve=500.)

all_sne=[sn2008ax,sn1996cb,sn1993J,sn2011dh,sn2011ei,sn2003bg]#,sn2011ei_without]

'''
for sn in all_sne:
    p.figure()
    x,y=sn.get_lightcurve()
    p.plot(x,log_n(y,mag_base),'bo')
    p.title(sn.sn_name)'''



for line in lines:
    line.interesting_sne=all_sne#this will plot data for all sne in all lines

def get_windows(band,sn,ts_days='none'):
    print '\n\n'
    print 'calculating window functions for %s in the %s band' %(sn.sn_name,band)
    windows=[[] for le in les]
    slit=Slit(0.,1.,subslits=8)
    if ts_days=='none':
        ts_days,_,_=sn.extract_data()#get days needed for window function
    sn.get_lightcurve(band)
    source=Source(sn.lightcurve[band],min(sn.lc_phases[band]),max(sn.lc_phases[band]),T-300*365.)
    for i in range(len(les)):
        print '\n'
        print 'Dust sheet %d of %d' %(i+1,len(les))
        les[i].fill_flux_array(source)
        window_function=les[i].window(slit)
        windows[i]=normalise(window_function(ts_days))
    print '\n'
    return windows


print '\n'


print 'producing spectra'


out={}
plot_data=[]
for line in lines:
    print 'Calculating for line %s' %line.name
    out[line.name]={'l0':line.l0}
    for sn in line.interesting_sne:
        print 'Supernova %s' %sn.sn_name
        ts_days,wlengths,_=sn.extract_data()#finds the range of phases and wavelthengths required
        sn.get_lightcurve()#this is used to see what band to use
        if line.band in sn.lc_flux.keys():
            band=line.band
        else:
            band=band_dictionary[line.band]
            print 'no %s band data. Using the %s band instead.' %(line.band,band)
            if band not in sn.lc_flux.keys():
                print 'Error! No %s band data. Ignoring this SN' %band
                continue
        
        full_lambdas=[]
        for a in wlengths:
            for b in a:
                full_lambdas.append(b)
        full_lambdas=np.sort(list(set(full_lambdas)))#remove duplicates and sort
        windows=get_windows(band,sn)
        windows.append(np.ones_like(windows[0]))#for unmodified
        lambdas=reduce_range(full_lambdas,line.minl,line.maxl)
        skip=skip_value(lambdas, lambda_bin)
        
        out[line.name][sn.sn_name]={'lambdas':lambdas[skip]}#for output
        
        spectra=[]
        mins=[]
        maxs=[]
        for i in range(len(windows)):
            spectra.append(apply_window(sn,windows[i], ts_days, lambdas,band=band))
            mins.append(apply_window(sn,windows[i], ts_days, lambdas,band=band,minmax='min'))
            maxs.append(apply_window(sn,windows[i], ts_days, lambdas,band=band,minmax='max'))
        
        labels=['le2116','le2521','le3923','unmodified']
        plot_spectra=[]
        if PLOTTING:
            p.figure()
            for i in range(len(spectra)):
                x,y,z=scale(spectra[i],mins[i],maxs[i])
                #x-=i*0.5
                #y-=i*0.5
                #z-=i*0.5
                p.plot(lambdas,x,colours[i],label=labels[i])
                #p.fill_between(lambdas,y,z,color=colours[i],alpha=0.4)
                p.plot(lambdas,y,colours[i]+'--')
                p.plot(lambdas,z,colours[i]+'--')
                out[line.name][sn.sn_name][labels[i]]=[x[skip],y[skip],z[skip]]
                plot_spectra.append([x,y,z])
            p.xlabel('Wavelength (A)')
            p.ylabel('Flux')
            p.title('line %s for %s' %(line_plotnames[line.name],sn.sn_name))
            p.legend()
        else:
            for i in range(len(spectra)):
                x,y,z=scale(spectra[i],mins[i],maxs[i])
                out[line.name][sn.sn_name][labels[i]]=[x[skip],y[skip],z[skip]]
                plot_spectra.append([x,y,z])
        plot_data.append([sn.sn_name,line.name,line.minl,line.maxl,lambdas,plot_spectra])

pdump(out,'individual.pkl')
pdump(plot_data,'plot_data/individual.pkl')
p.show()

