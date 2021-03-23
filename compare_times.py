'''
Created on Jul 11, 2015

@author: kieran
'''
from spectral_tools import flm_Template,apply_window,Pickle_Template
import numpy as np
from functions import log_n
from echo_tools_rect import *
import os
import pylab as p
from general_tools import progress_bar,pload,pdump
from copy import copy

T=2009.
d0=11000.#distance to Cas a (=3.4kpc
a=100.**(1./5)#base of the magnitude system
data_file='window_data/times_windows.dat'#file where the windows data is stored

colors=['b','g','r','c','m','k']#make it bigger if more ranges are added

times=[[-20,100],
       [-20,140],
       [-20,70]]

#                  sigma,   rho_0,      z0,         alpha,  distance
le_2116=dust_sheet(0.031,   627.9561,   437.11,     9.,     d0,     psf=0.97,delta=-39.97,offset=-0.42,profile=gaussian)#delta=-39.97,
le_2521=dust_sheet(0.093,   317.9613,   -9.89,      54.,    d0,     psf=1.13,delta=-0.91,offset=-0.01,profile=gaussian)#delta=-0.91,
le_3923=dust_sheet(0.614,   1203.8904,  2045.38,    7.,     d0,     psf=1.65,delta=-35.36,offset=-0.35,profile=gaussian)#delta=-35.36,


les=[le_2116,le_2521,le_3923]
le_labels=['le2116','le2521','le3923','unmodified']

band_dictionary={'r':'R_BVRI','i':'I_BVRI'}

#sn=flm_Template('sn1993J')
sn=Pickle_Template('meanspec/meanspec_GPinterp.p',continuum='one')
sn.sn_name='Template'
    
def scale(x):#scale such that the average value is 1
    x=np.array(x)
    norm=np.sum(x)/len(x)
    return x/norm

c=299.8#*10^3 km/s    
def lambda_to_v(l,l0):
    x=l0/l
    b=(1.-x)/(1.+x)
    return c*b


def get_label(rnge):
    minn,maxx=rnge
    return '%g<t<%g' %(minn,maxx)

def plot_compare(lambdas,data,labels,title):
    p.plot(lambdas,data[0],'k',lw=3,label='full_window')
    for i in range(1,len(data)):
        p.plot(lambdas,data[i],colors[i-1],label=get_label(labels[i-1]))
    p.legend()
    p.xlabel('Wavelength (A)')
    p.ylabel('Flux')
    p.title(title)
    
def plot_cross(xy,color):
    x,y=xy
    p.plot(x,y,'%sx'%color)
    
def find_min(x,y,rnge):
    minn,maxx=rnge
    i=0
    while x[i]<minn:
        i+=1
    lo=copy(i)
    while x[i]<maxx:
        i+=1
    hi=copy(i)
    a,b,c=np.polyfit(x[lo:hi],y[lo:hi],2)#fit is y=ax^2+bx+c
    xout=-b/(2*a)#has a minimum at x=-b/2a
    yout=a*xout**2+b*xout+c
    return (xout,yout)
    
def plot_mins(line,l0,ax):
    for key in line.keys():
        x,_=line[key]
        x=lambda_to_v(x,l0)
        ax.plot(0,x,'kx')
        ax.annotate(key,(0.1,x))
    ax.set_xlim(-1,3)
    
def get_windows(band,sn,recalculate=False):
    fname='window_data/%s_%s_band_window_data.dat' %(sn.sn_name,band)
    if recalculate:#forces recalculation
        os.remove(fname)        
    try:
        windows=pload(fname)
        print 'window functions read from file'
    except IOError:
        print '\n\n'
        print 'calculating window functions for %s in the %s band' %(sn.sn_name,band)
        windows=[[] for le in les]
        slit=Slit(0.,1.,subslits=8)
        ts_days,_,_=sn.extract_data()#get days needed for window function
        sn.get_lightcurve(band)
        source=Source(sn.lightcurve[band],min(sn.lc_phases[band]),max(sn.lc_phases[band]),T-300*365.)
        for i in range(len(les)):
            print '\n'
            print 'Dust sheet %d of %d' %(i+1,len(les))
            les[i].fill_flux_array(source)
            window_function=les[i].window(slit)
            windows[i]=normalise(window_function(ts_days))
        pdump(windows,fname)
        print '\n'
    return windows

class Line():
    def __init__(self,name,lambda_min,lambda_max,l0,band):
        self.name=name
        self.minl=lambda_min
        self.maxl=lambda_max
        self.band=band
        self.l0=l0
lines=[Line('H_alpha',5900.,6600.,6564.,'r'),
       Line('He',5300.,6000.,5876.,'V'),
       Line('Ca_III',7800.,8800.,8579.,'i'),
       Line('Fe_II',4900.,5300.,5169,'V')]


        
dle=1.
for line in lines:
    print 'Calculating for line %s' %line.name
    
    print 'calculating window functions'
    ts_days,wlengths,_=sn.extract_data()#finds the range of phases and wavelthengths required
    full_lambdas=[]
    for a in wlengths:
        for b in a:
            full_lambdas.append(b)
    full_lambdas=np.sort(list(set(full_lambdas)))#remove duplicates and sort
    lambdas=reduce_range(full_lambdas,line.minl,line.maxl)
    sn.get_lightcurve()#this is used to see what band to use
    if line.band in sn.lc_flux.keys():
        band=line.band
    else:
        band=band_dictionary[line.band]
        print 'no %s band data. Using the %s band instead.' %(line.band,band)
        if band not in sn.lc_flux.keys():
            print 'Error! No %s band data. Ignoring this SN' %band
            continue
    
    windows=get_windows(band,sn)
    test=np.array([len(windows[i])==len(ts_days) for i in range(len(windows))])
    if not test.all():
        print 'window functions out of date, need to recalculate'
        windows=get_windows(band,sn,recalculate=True)
    windows.append(np.ones_like(windows[0]))#for unmodified
    for i in range(len(windows)):
        full=apply_window(sn,windows[i], ts_days, lambdas,band=band)
        if line.name=='He' and i==0:
            p.plot(lambdas,scale(full)+i*dle,'k',lw=3,label='full_window')
        else:
            p.plot(lambdas,scale(full)+i*dle,'k',lw=3)
        for j in range(len(times)):
            rnge=times[j]
            tmin,tmax=rnge
            lo,hi=reduce_range(ts_days,tmin,tmax,indices=True)
            spectrum=apply_window(sn,windows[i][lo:hi], ts_days[lo:hi], lambdas,band=band)
            if line.name=='He' and i==0:
                p.plot(lambdas,scale(spectrum)+i*dle,colors[j],label=get_label(rnge))
            else:
                p.plot(lambdas,scale(spectrum)+i*dle,colors[j])
    
for line in lines:
    p.text((line.minl+line.maxl)/2.,4*dle+1,line.name)
for i in range(len(le_labels)):
    p.text(4500,i*dle+1,le_labels[i])
p.xlim(4480,max(full_lambdas))
p.legend()
p.show()