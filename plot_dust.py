'''
Created on Jun 30, 2015

creats plots in the rho-z plane of the dust sheet, ellipsoid and slit

@author: kieran
'''
import pylab as p
import numpy as np
from functions import *
from echo_tools import *
from general_tools import progress_bar
from copy import copy

'''use the following units
time: years
length: light years
angles: arcsecs (for the most part)
'''

c=1.

def plot_source(source,rhos,T):
    z=[]
    for dt in [source.left,0.,source.right]:
        t=T-(source.t0+dt)
        z.append(rhos**2/(2.*c*t)-(c*t)/2.)
    p.plot(z[1],rhos,'b')
    p.plot(z[0],rhos,'b--')
    p.plot(z[2],rhos,'b--')
    p.fill_betweenx(rhos, z[0], z[2],color='b',alpha=0.2)
    
def plot_dust(dust,rhos):
    z=[]
    rs=np.array([0.,dust.sigma/2.,dust.sigma])
    zds=rs/np.cos(dust.alpha)+dust.zd_0
    for zd in zds:
        z.append(zd-dust.a*rhos)
    p.plot(z[1],rhos,'g')
    p.plot(z[0],rhos,'g--')
    p.plot(z[2],rhos,'g--')
    p.fill_betweenx(rhos, z[0], z[2],color='g',alpha=0.2)
    
def plot_slit(slit,dust,zs):
    ones=np.ones_like(zs)
    minn=ones*dust.arcsec_to_ly(slit.rho_min+dust.offset)+dust.rho_0
    maxx=ones*dust.arcsec_to_ly(slit.rho_max+dust.offset)+dust.rho_0
    p.plot(zs,maxx,'r')
    p.plot(zs,minn,'r')
    p.fill_between(zs, minn,maxx,color='r',alpha=0.2)
    
def make_plot(dust,source,slit,rs,zs,T,title=False):
    rhos=copy(rs+dust.rho_0)
    plot_source(source, rhos, T)
    plot_dust(dust, rhos)
    plot_slit(slit,dust,zs)
    p.xlabel('z')
    p.ylabel('rho')
    p.xlim(min(zs),max(zs))
    p.ylim(min(rhos),max(rhos))
    if title:
        p.title(title)
    
    
    
    
    

T=0.#time of measurement. This is arbitrary, but fixes the origin of time.

dust_sheets=[]
source=Source(triangle(0.1,0.3),-0.1,0.3,T-300.)#paper says 300y old event 10000ly away but then wouldn't have seen it yet?
slit=Slit(0.,1.)



variables=[0.001,0.01,0.1,0.2,0.3,0.5]
rows=3
for sig in variables:
    dust_sheets.append(dust_sheet(sig,300.,0,45.,10000.,T))

'''
variables=[70.,45.,0.,-30.]
for alpha in variables:
    dust_sheets.append(dust_sheet(0.1,300.,0,alpha,10000.,T))
    
variables=[0,0.4,1.,2.]
sigmas=[0.001,0.1,0.3]
for sig in sigmas:
    for psf_sigma in variables:
        dust_sheets.append(dust_sheet(sig,300.,0,45.,10000.,T,psf=psf_sigma))

variables=[0.,-1.1,-2.2]
sigmas=[0.001,0.1,0.3]
for sig in sigmas:
    for ofst in variables:
        dust_sheets.append(dust_sheet(sig,300.,0,45.,10000.,T,offset=ofst))
'''


rho_arcs=np.linspace(-7,5,500)
rhos=dust_sheets[0].arcsec_to_ly(rho_arcs)#this assumes all dust sheets at the same distance. Change if necessary
zs=np.linspace(-1,1,2)

'''
for i in range(len(sigmas)):
    for j in range(len(variables)):
        p.figure()
        make_plot(dust_sheets[i*len(variables)+j],source,slit,rhos,zs,T,title='sigma=%g, offset=%g'%(sigmas[i],variables[j]))
'''

for i in range(len(variables)):
    p.figure()
    make_plot(dust_sheets[i],source,slit,rhos,zs,T,title='sigma=%g'%variables[i])


p.show()

