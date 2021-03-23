'''
Created on Jul 22, 2015

@author: kieran
'''
from spectral_tools import spec_plotter,Pickle_Template,apply_window,get_band
from functions import log_n
from echo_tools_rect import *
import pylab as p
import os
from general_tools import progress_bar,pload,pdump,reduce_range
import sys
import matplotlib.lines as mlines

import matplotlib
matplotlib.rcParams.update({'font.size': 16})#increase font size

T=2009.*365
d0=11000.#distance to Cas a (=3.4kpc
mag_base=100.**(1./5)#base of the magnitude system
if 'data' in sys.argv:
    PLOTTING=False
else:
    PLOTTING=True


def scale(x,er1,er2):#scale such that the average value is 1, er1 and er2 are errors that are scaled the same way
    x=np.array(x)
    norm=np.abs(np.sum(x)/len(x))
    return (x/norm,er1/norm,er2/norm)

class Line():
    def __init__(self,name,lambda_min,lambda_max,l0,band):
        self.name=name
        self.minl=lambda_min
        self.maxl=lambda_max
        self.band=band
        self.l0=l0

#                  sigma,   rho_0,      z0,         alpha,  distance
le_2116=dust_sheet(0.031,   627.9561,   437.11,     9.,     d0,     psf=0.97,delta=-39.97,offset=-0.42,profile=gaussian)#delta=-39.97,
le_2521=dust_sheet(0.093,   317.9613,   -9.89,      54.,    d0,     psf=1.13,delta=-0.91,offset=-0.01,profile=gaussian)#delta=-0.91,
le_3923=dust_sheet(0.614,   1203.8904,  2045.38,    7.,     d0,     psf=1.65,delta=-35.36,offset=-0.35,profile=gaussian)#delta=-35.36,

les=[le_2116,le_2521,le_3923]
colours=['r','b','y','k']#last colour is for unmodified

lines=[Line('H_alpha',5900.,6600.,6564.,'r'),
       Line('He',5350.,6000.,5876.,'V'),
       Line('He',5350.,6000.,5876.,'r'),
       Line('Ca_III',7800.,8800.,8579.,'i')]

       #Line('Fe_II',4850.,5350.,5169,'V')]]


bands=['B','V','r','i']
line_colors={'H_alpha':'b','He':'c','Fe_II':'g','Ca_III':'r'}
band_colors=['b','g','r','m']
line_plotnames={'He':'He 5876','H_alpha':r'$H\alpha$','Ca_III':'Ca II','Fe_II':'Fe 5169'}
labels=['le2116','le2521','le3923','unmodified']

x=np.linspace(3000,10000,1000)
if PLOTTING:
    plot_bands=[]
    for i in range(len(bands)):
        band=bands[i]
        y=get_band(band)(x)
        p.plot(x,y,label=band,color=band_colors[i])
        plot_bands.append([band,x,y])
    line_names=[]
    plot_lines=[]
    for line in lines:
        if line.name not in line_names:
            p.axvline(line.minl,color=line_colors[line.name])
            p.axvline(line.maxl,color=line_colors[line.name])
            p.axvspan(line.minl,line.maxl,alpha=0.4,label=line_plotnames[line.name],color=line_colors[line.name])
            line_names.append(line.name)
            plot_lines.append([line.name,line.minl,line.maxl])
    p.legend(loc='lower left')
    p.xlabel('Wavelength (A)')
    p.ylabel('Flux')
    pdump({'bands':plot_bands,'lines':plot_lines},'plot_data/bands_lines.pkl')


sn=Pickle_Template('meanspec/meanspec_GPinterp.p',continuum='one')
sn.get_lightcurve()
ts_days,wlengths,_=sn.extract_data()#finds the range of phases and wavelengths required

full_lambdas=[]
for a in wlengths:
    for b in a:
        full_lambdas.append(b)
full_lambdas=np.sort(list(set(full_lambdas)))#remove duplicates and sort

def get_windows(band,recalculate=False):
    fname='window_data/%s_band_window_data.dat' %band
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


print '\n'


print 'producing spectra'

out={}
plot_data=[]
for line in lines:
    windows=get_windows(line.band)
    test=np.array([len(windows[i])==len(ts_days) for i in range(len(windows))])
    if not test.all():
        print 'window functions out of date, need to recalculate'
        windows=get_windows(line.band,recalculate=True)
    windows.append(np.ones_like(windows[0]))#for unmodified
    lambdas=reduce_range(full_lambdas,line.minl,line.maxl)
    
    out[line.name]={'lambdas':lambdas,'l0':line.l0}#for output
    
    spectra=[]
    mins=[]
    maxs=[]
    for i in range(len(windows)):
        spectra.append(apply_window(sn,windows[i], ts_days, lambdas,band=line.band))
        mins.append(apply_window(sn,windows[i], ts_days, lambdas,band=line.band,minmax='min'))
        maxs.append(apply_window(sn,windows[i], ts_days, lambdas,band=line.band,minmax='max'))
    
    
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
            out[line.name][labels[i]]=[x,y,z]
            plot_spectra.append([x,y,z])
        #x,y=find_min(lambdas,spectra[-1])
        #p.axvline(x,color='k',ls='--')
            
        #p.legend(loc='lower right')
        if line.name=='Ca_III':
            p.legend(loc='upper left')
            p.ylim(0.55,1.6)
        elif line.name=='H_alpha':
            a=mlines.Line2D([],[],color='k',label='Mean')
            b=mlines.Line2D([],[],color='k',ls='--',label=r'$\pm1\sigma$')
            p.legend(handles=[a,b],loc='lower right')
            p.ylim(0.7,1.35)
        p.xlim(line.minl-50,line.maxl+50)
        ax=p.gca()
        if line.name=='H_alpha':
            y_text=0.9
        else:
            y_text=0.05
        ax.text(.5,y_text,'%s line, calibrated with %s band' %(line_plotnames[line.name],line.band),horizontalalignment='center',transform=ax.transAxes)
        p.xlabel('Wavelength (A)')
        p.ylabel('Flux')
        p.savefig('spectra_%s_%s' %(line.name,line.band))
    else:
        for i in range(len(spectra)):
            x,y,z=scale(spectra[i],mins[i],maxs[i])
            out[line.name][labels[i]]=[x,y,z]
            plot_spectra.append([x,y,z])
    plot_data.append([line.name,line.band,line.minl,line.maxl,lambdas,plot_spectra])
pdump(out,'features_le.pkl')
pdump(plot_data,'plot_data/features.pkl')
p.show()




