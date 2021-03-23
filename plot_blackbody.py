'''
Created on Jul 13, 2015

@author: kieran
'''
import pylab as p
import numpy as np


kb=8.6173324/(1e5)#eV/K
hbarc=1973.27#eVA
T=7000.#K

def Planck(e):
    pre=e**5
    post=np.exp(1./(e*kb*T))-1
    return 1./(pre*post)

x=np.linspace(0,15000,1000)
Es=x/hbarc
y=Planck(Es)
y=y/np.nanmax(y)
p.plot(x,y)
p.yscale('log')
p.ylim(1./1e5,10)
p.show()