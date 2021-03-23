'''
Created on Jul 3, 2015

@author: kieran
'''
from spectral_tools import Sav_Template,apply_window
from functions import log_n
from echo_tools import *
import pylab as p
from general_tools import progress_bar,pload,pdump
from time import time
import pickle as pkl

T=2009.
d0=11000.#distance to Cas a (=3.4kpc
mag_base=100.**(1./5)#base of the magnitude system

WINDOW_FNAME='template_windows.dat'

def scale(x):#scale such that the average value is 1
    x=np.array(x)
    norm=np.abs(np.sum(x)/len(x))
    return x/norm

#                  sigma,   rho_0,      z0,         alpha,  distance
le_2116=dust_sheet(0.031,   627.9561,   437.11,     9.,     d0,     T,psf=0.97,offset=-0.42,delta=-39.97,profile=gaussian)
le_2521=dust_sheet(0.093,   317.9613,   -9.89,      54.,    d0,     T,psf=1.13,offset=-0.01,delta=-0.91,profile=gaussian)
le_3923=dust_sheet(0.614,   1203.8904,  2045.38,    7.,     d0,     T,psf=1.65,offset=-0.35,delta=-35.36,profile=gaussian)

les=[le_2116,le_2521,le_3923]
colours=['r','b','y']

sn=Sav_Template('IIb_1specperSN')
sn.get_lightcurve()
ts_days,_,_=sn.extract_data()#finds the range of phases required
#to add ts_days=np.array(x,y,500)
sn.continuum='cubic'#adds a cubic continuum
sn.get_lightcurve()
mags=sn.lightcurve['V'](ts_days)
mags=normalise_mags(mags)


try:
    windows=pload(WINDOW_FNAME)
    print 'window functions read from file'
except IOError:
    print '\n\n'
    print 'calculating window functions'
    windows=[[] for le in les]
    slit=Slit(0.,1.,subslits=8)
    ts_years=ts_days/365.#todo change echo_tools to work in days
    source=Source(sn.lightcurve['V'],min(sn.lc_phases)/365.,max(sn.lc_phases/365.),T-300)
    for i in range(len(les)):
        print '\n'
        print 'Dust sheet %d of %d' %(i+1,len(les))
        for j in range(len(ts_years)):
            progress_bar(j,len(ts_years))
            windows[i].append(les[i].window(ts_years[j],T,slit,source=source))
        windows[i]=normalise(np.array(windows[i]))
    pdump(windows,WINDOW_FNAME)



        
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
ax[0].plot(ts_days,np.ones_like(ts_days),'k')
for i in range(len(windows)):
    ax[0].plot(ts_days,normalise(windows[i]),colours[i])
#ax[0].plot([-10,-10],[0,1.2],'k--')
#ax[0].plot([40,40],[0,1.2],'k--')
ax[0].set_ylim(0,1.2)
ax[0].set_ylabel('Window Function')


ax[1].plot(ts_days,normalise_mags(mags),'k')
for i in range(len(light_curves)):
    ax[1].plot(ts_days,light_curves[i],colours[i])
#ax[1].plot([-10,-10],[6,-1],'k--')
#ax[1].plot([40,40],[6,-1],'k--')
ax[1].set_ylim(6,-1)

ax[1].set_ylabel('V Magnitude')

ax[1].set_xlim(min(ts_days),250)
ax[1].set_xlabel('Time, days')
ax[0].legend(['sn1993j','le2116','le2521','le3923'])




p.figure()
#p.plot(lambdas,remove_continuum(lambdas,unmodified),'k')
p.plot(lambdas,scale(unmodified),'k')
for i in range(len(les)):
    #p.plot(lambdas,remove_continuum(lambdas,spectra[i]),colours[i])
    p.plot(lambdas,scale(spectra[i]),colours[i])

wfspectra = np.array(lambdas, scale(spectra[0]), scale(spectra[1]), scale(spectra[2]))
f = open ('Templates_lewf.pkl','wb')
pl.dump(wfspectra,f)
p.legend(['unmodified','le2116','le2521','le3923'],loc='lower right')
p.xlabel('Wavelength (A)')
p.ylabel('Flux')
p.title('Type IIb mean spectrum')
p.savefig('%s.png'%sn.name)


p.show()




