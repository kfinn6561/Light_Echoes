'''
Created on Jun 17, 2015

@author: kieran
'''
from spectral_tools import spec_plotter,apply_window
from functions import log_n
from echo_tools import *
import pylab as p
from general_tools import progress_bar,pload,pdump

T=2009.
d0=11000.#distance to Cas a (=3.4kpc
a=100.**(1./5)#base of the magnitude system

def scale(x):#scale such that the average value is 1
    x=np.array(x)
    norm=np.sum(x)/len(x)
    return x/norm

#                  sigma,   rho_0,      z0,         alpha,  distance
le_2116=dust_sheet(0.031,   627.9561,   437.11,     9.,     d0,     T,psf=0.97,offset=-0.42,delta=-39.97,profile=gaussian)
le_2521=dust_sheet(0.093,   317.9613,   -9.89,      54.,    d0,     T,psf=1.13,offset=-0.01,delta=-0.91,profile=gaussian)
le_3923=dust_sheet(0.614,   1203.8904,  2045.38,    7.,     d0,     T,psf=1.65,offset=-0.35,delta=-35.36,profile=gaussian)

les=[le_2116,le_2521,le_3923]
colours=['r','b','y']

  
print('Loading the data')
sn1993j=spec_plotter()
sn1993j.loadspeclist('sn1993j')
sn1993j.specreductiontype='dered-warp'
sn1993j.loadspectra()
sn1993j.loadlc_MLCS()

ts,mags=sn1993j.get_lightcurve()#this will need to be updated when there are more than one SN

ts=ts/365#convert to years
mags=normalise_mags(mags)


'''
window_left=-50.#days
window_right=250.#days
ts=np.linspace(window_left/365.,window_right/365.,500)
'''


try:
    windows=pload('rest_windows.dat')
    print 'window functions read from file'
except IOError:
    print '\n\n'
    print 'calculating window functions'
    slit=Slit(0.,1.,subslits=8)
    #slit=Slit(0.,1.)
    source=Source(sn1993j.lightcurve,min(sn1993j.lc_phases)/365.,max(sn1993j.lc_phases/365.),T-300)  
    windows=[[] for le in les]
    for i in range(len(les)):
        print '\n'
        print 'Dust sheet %d of %d' %(i+1,len(les))
        for j in range(len(ts)):
            progress_bar(j,len(ts))
            windows[i].append(les[i].window(ts[j],T,slit,source=source))
        windows[i]=normalise(np.array(windows[i]))
    pdump(windows,'rest_windows.dat')


        
print '\n'
print 'convolving with the light curves'
light_curves=[]
for i in range(len(les)):
    light_curves.append(mags-log_n(windows[i],a))


f,ax=p.subplots(2,1,sharex=True)
ax[0].plot(ts*365,np.ones_like(ts),'k')
for i in range(len(windows)):
    ax[0].plot(ts*365,normalise(windows[i]),colours[i])
#ax[0].plot([-10,-10],[0,1.2],'k--')
#ax[0].plot([40,40],[0,1.2],'k--')
ax[0].set_ylim(0,1.2)
ax[0].set_ylabel('Window Function')


ax[1].plot(ts*365,normalise_mags(mags),'k')
for i in range(len(light_curves)):
    ax[1].plot(ts*365,light_curves[i],colours[i])
#ax[1].plot([-10,-10],[6,-1],'k--')
#ax[1].plot([40,40],[6,-1],'k--')
ax[1].set_ylim(6,-1)

ax[1].set_ylabel('Magnitude')

ax[1].set_xlim(min(ts)*365,250)
ax[1].set_xlabel('Time')
ax[0].legend(['sn1993j','le2116','le2521','le3923'])


print 'producing spectra'
phase,wlength,flux=sn1993j.extract_data()

lambdas=np.linspace(4800,9000,500)
spectra=[]
for i in range(len(les)):
    spectra.append(apply_window(sn1993j,windows[i], ts*365, lambdas))
unmodified=apply_window(sn1993j,np.ones_like(windows[0]), ts*365, lambdas)#unmodified spectrum is weighted with flat window fn


p.figure()
p.plot(lambdas,scale(unmodified),'k')#may need to apply some normalisation
for i in range(len(les)):
    p.plot(lambdas,scale(spectra[i]),colours[i])
    
p.legend(['sn1993j','le2116','le2521','le3923'])
p.xlabel('Wavelength (A)')
p.ylabel('Flux')

p.show()






