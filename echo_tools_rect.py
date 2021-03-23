'''
Created on June 2nd, 2015

RUN SOURCE EXPORT.SOURCEME BEFORE STARTING


@author: kieran
'''
import numpy as np
from functions import *
from copy import copy
from general_tools import find_nearest,reduce_range
from scipy.interpolate import interp1d


'''use the following units
time: days
length: light years
angles: arcsecs (for the most part)
'''


c=1./365. #ly/d speed of light
Nsigma=4.
rhos=np.linspace(-7,5,201)
tres=301#accuracy of t integrals

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
    def __init__(self,sigma,rho_0,z0,alpha,distance,psf=0,offset=0,delta=0,profile=get_tophat):
        self.sigma=sigma
        self.rho_0=rho_0
        self.z0=z0
        self.alpha=np.deg2rad(alpha)#alpha is in degrees
        self.D=profile(sigma)
        self.a=np.tan(self.alpha)
        self.zd_0=self.z0+self.a*self.rho_0
        self.distance=np.sqrt((distance-z0)**2+rho_0**2)
        if psf!=0:
            self.psf_sigma=psf/(2.*np.sqrt(2.*np.log(2.)))#covert FWHM to variance
            self.psf=gaussian(self.psf_sigma)
        else:
            self.psf_sigma=0
        self.offset=offset
        self.delta=np.deg2rad(delta)
        '''Not using T0, should still work'''
        self.rhos=rhos
        self.drho=rhos[1]-rhos[0] 
        self.rhos_ly=self.arcsec_to_ly(self.rhos)
        self.t0=(-self.z0+np.sqrt(self.z0**2+self.rho_0**2))/c#time delay of centre of explosion
        
        
    
    def get_r(self,rho,t):#rho in ly
        zd=(rho**2)/(2*c*t)-(c*t)/2.+self.a*rho
        return np.cos(self.alpha)*(zd-self.zd_0)
    
    def flux(self,rho,t):#rho in ly
        r=self.get_r(rho+self.rho_0, self.t0-t)
        return self.D(r)*self.source(t)
    
    def fill_flux_array(self,source):
        self.source=source
        self.ts=np.linspace(source.left,source.right,tres)
        self.dt=self.ts[1]-self.ts[0]        
        self.flux_array=np.ones([len(self.ts),len(self.rhos)])*self.dt*self.drho
        
        self.t_array=np.array([self.ts for i in range(len(self.rhos))]).T
        self.rho_ly_array=np.array([self.rhos_ly for i in range(len(self.ts))])
        
        out_flux=self.flux(self.rho_ly_array,self.t_array)
        self.flux_array*=out_flux
        self.apply_seeing()
    
    
    def apply_seeing(self):
        if self.psf_sigma==0:
            return
        out=[]
        for i in range(len(self.rhos)):
            out.append(np.outer(self.flux_array.T[i],self.psf(self.rhos-self.rhos[i])))
        self.flux_array*=(1-self.psf(0))#cancels out the above effect applying psf to same square
        for m in out:
            self.flux_array+=m
    
    def profile(self):
        try:
            out= self.profile_vals
        except AttributeError:
            self.profile_vals=np.sum(self.flux_array,0)
            out= self.profile_vals
        return interp1d(self.rhos,out,bounds_error=False,fill_value=0)
    
    def total_flux(self,slit):
        rho_min=slit.rho_min+self.offset
        rho_max=slit.rho_max+self.offset
        lo,hi=reduce_range(self.rhos,rho_min,rho_max,indices=True)
        return np.sum(self.flux_array[:,lo:hi])
        
    def effective_lcv(self,slit):
        rho_min=slit.rho_min+self.offset
        rho_max=slit.rho_max+self.offset
        lo,hi=reduce_range(self.rhos,rho_min,rho_max,indices=True)
        return np.sum(self.flux_array[:,lo:hi],1)
    
    def window(self,slit):
        slit.sigma/=np.cos(self.delta)
        slit.subsig/=np.cos(self.delta)
        if len(slit.us)==1:
            out=self.subwindow(slit)
        else:
            weights=[]
            windows=[]
            for u in slit.us:
                subslit=Slit(u*np.sin(self.delta),slit.subsig)
                windows.append(self.subwindow(subslit))
                weights.append(self.total_flux(subslit))
            weights=np.array(weights)
            weights=np.outer(weights,np.ones_like(windows[0]))
            windows=np.array(windows)
            out=np.sum(weights*windows,0)/np.sum(weights,0)
        return interp1d(self.ts,out,bounds_error=False,fill_value=0)
    
    def subwindow(self,slit):
        return self.effective_lcv(slit)/self.source(self.ts)
    
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





