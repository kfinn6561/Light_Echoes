'''
Created on Jul 3, 2015

@author: kieran
'''
from spectral_tools import Pickle_Template,apply_window
from functions import log_n
from echo_tools_rect import *
import pylab as p
from general_tools import progress_bar,pload,pdump,increase_sampling
from time import time

import matplotlib
matplotlib.rcParams.update({'font.size': 16})#increase font size


T=2009.*365.
d0=11000.#distance to Cas a (=3.4kpc
mag_base=100.**(1./5)#base of the magnitude system

WINDOW_FNAME='window_data/continuum_template_windows.dat'

def scale(x):#scale such that the average value is 1
    x=np.array(x)
    norm=np.abs(np.sum(x)/len(x))
    return x/norm

#                  sigma,   rho_0,      z0,         alpha,  distance
le_2116=dust_sheet(0.031,   627.9561,   437.11,     9.,     d0,     psf=0.97,delta=-39.97,offset=-0.42,profile=gaussian)#delta=-39.97,
le_2521=dust_sheet(0.093,   317.9613,   -9.89,      54.,    d0,     psf=1.13,delta=-0.91,offset=-0.01,profile=gaussian)#delta=-0.91,
le_3923=dust_sheet(0.614,   1203.8904,  2045.38,    7.,     d0,     psf=1.65,delta=-35.36,offset=-0.35,profile=gaussian)#delta=-35.36,

les=[le_2116,le_2521,le_3923]
colours=['r','b','y']

sn=Pickle_Template('meanspec/meanspec_GPinterp.p')
sn.get_lightcurve()
ts_days,_,_=sn.extract_data()#finds the range of phases required
#ts_days=increase_sampling(ts_days,10)
sn.get_lightcurve()
mags=sn.lightcurve_mags['V'](ts_days)
mags=normalise_mags(mags)


print '\n\n'
print 'calculating window functions'
windows=[[] for le in les]
slit=Slit(0.,1.,subslits=8)
source=Source(sn.lightcurve['V'],min(sn.lc_phases['V']),max(sn.lc_phases['V']),T-300*365.)
for i in range(len(les)):
    print '\n'
    print 'Dust sheet %d of %d' %(i+1,len(les))
    les[i].fill_flux_array(source)
    window_function=les[i].window(slit)
    windows[i]=normalise(window_function(ts_days))



        
print '\n'


print 'producing spectra'


phase,wlength,flux=sn.extract_data()

lambdas=np.linspace(4200,8500,500)
spectra=[]
for i in range(len(les)):
    spectra.append(apply_window(sn,windows[i], ts_days, lambdas))
unmodified=apply_window(sn,np.ones_like(windows[0]), ts_days, lambdas)#unmodified spectrum is weighted with flat window fn


print '\n'
print 'convolving with the light curves'
light_curves=[]
for i in range(len(les)):
    light_curves.append(mags-log_n(windows[i],mag_base))


f,ax=p.subplots(2,1,sharex=True)
p.tight_layout(h_pad=-1)
wfs=[np.ones_like(ts_days)]
ax[0].plot(ts_days,np.ones_like(ts_days),'k')
for i in range(len(windows)):
    y=normalise(windows[i])
    ax[0].plot(ts_days,y,colours[i])
    wfs.append(y)
#ax[0].plot([-10,-10],[0,1.2],'k--')
#ax[0].plot([40,40],[0,1.2],'k--')
ax[0].set_ylim(0,1.2)
ax[0].set_ylabel('Window Function')


lcs=[normalise_mags(mags)]
ax[1].plot(ts_days,normalise_mags(mags),'k')
for i in range(len(light_curves)):
    ax[1].plot(ts_days,light_curves[i],colours[i])
    lcs.append(light_curves[i])
#ax[1].plot([-10,-10],[6,-1],'k--')
#ax[1].plot([40,40],[6,-1],'k--')
ax[1].set_ylim(6,-1)
ax[1].set_yticks([6,5,4,3,2,1,0])

ax[1].set_ylabel('V Magnitude')

ax[1].set_xlim(min(ts_days),max(ts_days))
ax[1].set_xlabel('Time (days)')
ax[1].legend(['unmodified','le2116','le2521','le3923'],loc='lower left')

plot_out={'windows':[ts_days,wfs],'light_curves':[ts_days,lcs]}
pdump(plot_out,'plot_data/window_functions.pkl')


p.figure()
#p.plot(lambdas,remove_continuum(lambdas,unmodified),'k')
p.plot(lambdas,scale(unmodified),'k')
for i in range(len(les)):
    #p.plot(lambdas,remove_continuum(lambdas,spectra[i]),colours[i])
    p.plot(lambdas,scale(spectra[i]),colours[i])
    
p.legend(['unmodified','le2116','le2521','le3923'],loc='upper right')
p.xlabel('Wavelength (A)')
p.ylabel('Flux')
p.title('Type IIb mean spectrum')
p.savefig('continuum_template.png')


p.show()




