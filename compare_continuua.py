'''
Created on Aug 19, 2015

@author: kieran
'''
from spectral_tools import Pickle_Template,apply_window
from functions import log_n
from echo_tools_rect import *
import pylab as p
from general_tools import progress_bar,pload,pdump,increase_sampling,reduce_range
import os
from copy import deepcopy

import matplotlib
matplotlib.rcParams.update({'font.size': 16})#increase font size

T=2009.
d0=11000.#distance to Cas a (=3.4kpc
mag_base=100.**(1./5)#base of the magnitude system

def scale(x,er1,er2):#scale such that the average value is 1, er1 and er2 are errors that are scaled the same way
    x=np.array(x)
    norm=np.abs(np.sum(x)/len(x))
    return (x/norm,er1/norm,er2/norm)

#                  sigma,   rho_0,      z0,         alpha,  distance
le_2116=dust_sheet(0.031,   627.9561,   437.11,     9.,     d0,     psf=0.97,offset=-0.42,delta=-39.97,profile=gaussian)
le_2521=dust_sheet(0.093,   317.9613,   -9.89,      54.,    d0,     psf=1.13,offset=-0.01,delta=-0.91,profile=gaussian)
le_3923=dust_sheet(0.614,   1203.8904,  2045.38,    7.,     d0,     psf=1.65,offset=-0.35,delta=-35.36,profile=gaussian)

les=[le_2116,le_2521,le_3923]
le_names=['le2116','le2521','le3923']

lines=pload('features_le.pkl')
line_names=['He','H_alpha','Ca_III']
cols=len(line_names)

sn_conti=Pickle_Template('meanspec/meanspec_GPinterp_conti.p')
sn_conti.plotname='Measured Continuum'
sn_cubic=Pickle_Template('meanspec/meanspec_GPinterp.p',continuum='cubic')
sn_cubic.plotname='Cubic Continuum'
sn_flat=Pickle_Template('meanspec/meanspec_GPinterp.p',continuum='one')

sne=[sn_conti,sn_cubic]
for band in ['V','r','i']:
    SN=deepcopy(sn_flat)
    SN.plotname='%s band' %band
    SN.band=band
    sne.append(SN)


ts_days,_,_=sn_conti.extract_data()#finds the range of phases required, assumes same for all sne


def get_windows(band,recalculate=False):
    fname='window_data/%s_band_compare.dat' %band
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
        sn_conti.get_lightcurve(band)
        source=Source(sn_conti.lightcurve[band],min(sn_conti.lc_phases[band]),max(sn_conti.lc_phases[band]),T-300)
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

lambdas=np.linspace(4200,9000,1000)
f,ax=p.subplots(len(les),1,sharex=True)
color_cycle = ax[0]._get_lines.color_cycle
spectra={sn.plotname:{} for sn in sne}
for sn in sne:
    color=next(color_cycle)
    sn.color=color
    try:
        band=sn.band
    except AttributeError:
        band='V'
        
    phase,wlength,flux=sn.extract_data()
    sn.get_lightcurve()
    mags=sn.lightcurve[band](ts_days)
    mags=normalise_mags(mags)
    windows=get_windows(band)
    test=np.array([len(windows[i])==len(ts_days) for i in range(len(windows))])
    if not test.all():
        print 'window functions out of date, need to recalculate'
        windows=get_windows(band,recalculate=True)

    for i in range(len(les)):
        spectrum=apply_window(sn,windows[i], ts_days, lambdas,band=band)
        spectra[sn.plotname][le_names[i]]=spectrum
        minn=apply_window(sn,windows[i], ts_days, lambdas,band=band,minmax='min')
        maxx=apply_window(sn,windows[i], ts_days, lambdas,band=band,minmax='max')
        #spectrum,minn,maxx=scale(spectrum,minn,maxx)
        ax[i].plot(lambdas,spectrum,color=color,label=sn.plotname)
        #ax[i].plot(lambdas,minn,color=color,ls='--')
        #ax[i].plot(lambdas,maxx,color=color,ls='--')
for i in range(len(les)):
    ax[i].set_ylabel(le_names[i])
    ax[i].get_yaxis().set_ticks([])
ax[-1].set_xlabel('Wavelength (A)')
ax[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=len(sne))
b,t=ax[0].get_ylim()
w=t-b
ax[0].set_ylim(b,b+1.3*w)#increase the width so that the legend can be seen

f_split,ax_split=p.subplots(len(lines),len(les),sharex='col',)
handles=[]
labels=[]
for i in range(len(line_names)):
    for j in range(len(le_names)):
        lamrange=lines[line_names[i]]['lambdas']
        lo,hi=reduce_range(lambdas,lamrange[0],lamrange[-1],indices=True)
        x=lambdas[lo:hi]
        for sn in sne:#V and i
            y=spectra[sn.plotname][le_names[j]][lo:hi]
            y,_,_=scale(y,1,1)
            temp=ax_split[j][i].plot(x,y,color=sn.color,label=sn.plotname)
            if i==0 and j==0:
                handles.append(temp[0])
                labels.append(sn.plotname)
            ax_split[j][i].get_yaxis().set_ticks([])
        
for i in range(len(line_names)):
    ax_split[0][i].set_title(line_names[i])
    b,t=ax_split[0][i].get_ylim()
    w=t-b
    ax_split[0][i].set_ylim(b,b+1.3*w)#increase the width so that the legend can be seen
for j in range(len(le_names)):
    ax_split[j][0].set_ylabel(le_names[j])
f_split.legend(handles,labels,loc='upper center', ncol=len(sne))




p.show()




