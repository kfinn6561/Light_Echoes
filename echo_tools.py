'''
Created on June 2nd, 2015

RUN SOURCE EXPORT.SOURCEME BEFORE STARTING


@author: kieran
'''
import numpy as np
from functions import *
from copy import copy

'''use the following units
time: days
length: light years
angles: arcsecs (for the most part)
'''

c=1./365. #ly/d speed of light
Nsigma=4.

class Source():
    def __init__(self,f,left,right,t0):
        self.f=f
        self.left=left
        self.right=right
        self.t0=t0
    def __call__(self,x):
        return self.f(x)
    
class Slit():
    def __init__(self,offset,sigma,subslits=1):
        self.offset=offset
        self.sigma=sigma
        self.rho_min=offset-sigma/2.
        self.rho_max=offset+sigma/2.
        us=np.linspace(-sigma/2,sigma/2,subslits+1)
        self.us=(us[:-1]+us[1:])/2
        self.subsig=sigma/subslits #I think this is actually wrong
        #self.subsig=sigma

def get_tophat(sigma,centre=0.):
    return lambda x:top_hat(x,centre-sigma/2.,centre+sigma/2.)


class dust_sheet():
    def __init__(self,sigma,rho_0,z0,alpha,distance,T0,psf=0,offset=0,delta=0,profile=get_tophat):
        self.sigma=sigma
        self.rho_0=rho_0
        self.z0=z0
        self.alpha=np.deg2rad(alpha)#alpha is in degrees
        self.D=profile(sigma)
        self.a=np.tan(self.alpha)
        self.zd_0=self.z0+self.a*self.rho_0
        self.distance=np.sqrt((distance-z0)**2+rho_0**2)
        if psf!=0:
            self.psf_sigma=self.arcsec_to_ly(psf)#psf is in arcsecs
            self.psf_sigma=self.psf_sigma/(2.*np.sqrt(2.*np.log(2.)))#covert FWHM to variance
            self.psf=gaussian(self.psf_sigma)
        else:
            self.psf_sigma=0
        self.offset=offset
        self.delta=np.deg2rad(delta)
        dt0=(-self.z0+np.sqrt(self.z0**2+self.rho_0**2))/c#time delay of centre of explosion
        self.t0=T0-dt0#time of peak of explosion, expect T0 in years
    
    def get_r(self,rho,t):#rho in ly
        zd=(rho**2)/(2*c*t)-(c*t)/2.+self.a*rho
        return np.cos(self.alpha)*(zd-self.zd_0)
    
    def flux(self,rho,t,T,source):#rho in ly
        r=self.get_r(rho, T-t)
        return self.D(r)*source(t-self.t0)
    
    def profile(self,rho,T,source):
        if self.psf_sigma==0:
            out=integrate(lambda x:self.flux(rho,x+self.t0,T,source),source.left,source.right,N=2000)
        else:
            f=lambda x,dt:self.flux(rho+x,dt+self.t0,T,source)
            out=Integrate(f,[-Nsigma*self.psf_sigma, source.left], [Nsigma*self.psf_sigma, source.right],method='cuhre')
        return out
    
    def total_flux(self,T,source,slit):
        rho_min=self.arcsec_to_ly(slit.rho_min+self.offset)
        rho_max=self.arcsec_to_ly(slit.rho_max+self.offset)
        if self.psf_sigma==0:
            f=lambda x,y:self.flux(x+self.rho_0,y+self.t0,T,source)
            return Integrate(f,[rho_min,source.left],[rho_max,source.right],method='cuhre')#may want to change this
        else:
            f=lambda x,y,z:self.flux(x+z+self.rho_0,y+self.t0,T,source)*self.psf(z)
            a=[rho_min,source.left,-Nsigma*self.psf_sigma]
            b=[rho_max,source.right,Nsigma*self.psf_sigma]
            return Integrate(f,a,b,method='cuhre')#3D integration
    
    def window(self,dt,T,slit,source=None):
        slit.sigma/=np.cos(self.delta)
        slit.subsig/=np.cos(self.delta)
        if len(slit.us)==1:
            return self.subwindow(dt,T,slit)
        else:
            weights=[]
            windows=[]
            for u in slit.us:
                subslit=Slit(u*np.sin(self.delta),slit.subsig)
                windows.append(self.subwindow(dt,T,subslit))
                weights.append(self.total_flux(T,source,subslit))
            weights=np.array(weights)
            windows=np.array(windows)
            out=np.sum(weights*windows)/np.sum(weights)
            return out
    
    def subwindow(self,dt,T,slit):
        left=self.arcsec_to_ly(slit.rho_min+self.offset)
        right=self.arcsec_to_ly(slit.rho_max+self.offset)
        t=dt+self.t0
        if self.psf_sigma==0:
            out=integrate(lambda x: self.D(self.get_r(x+self.rho_0,T-t)),left,right,N=1000)
        else:
            f=lambda x,y:self.D(self.get_r(x+self.rho_0+y,T-t))*self.psf(y)
            #out=integrate_2d(f,-3*self.psf_sigma,3*self.psf_sigma,left,right,Nx=1000,Ny=100)
            out=Integrate(f,[-Nsigma*self.psf_sigma,left],[Nsigma*self.psf_sigma,right],method='cuhre')
        return out
    
    def ly_to_arcsec(self,x):
        theta=x/(self.distance-self.z0)#radians
        theta=np.rad2deg(theta)#degrees
        theta=theta*3600.#arcsecs
        return theta
    def arcsec_to_ly(self,theta):
        theta=np.deg2rad(theta/3600.)
        x=theta*(self.distance-self.z0)
        return x
               
        
def normalise(x):#normalises so the maximum is 1
    x=np.array(x)
    return x/np.nanmax(x)#nanmax ignores nans when calculating the max

def normalise_mags(x):#normalises magnitudes so that the minimum is 0
    x=np.array(x)
    return x-np.nanmin(x)





