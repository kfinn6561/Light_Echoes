'''
Created on Jul 16, 2015

@author: kieran
'''
from scipy.interpolate import interp1d,InterpolatedUnivariateSpline,UnivariateSpline
import numpy as np
import pylab as p
from copy import copy
from functions import fit_spline,plot_spline

def test_splines(x,y):
    p.figure()
    p.plot(x,y,'ro')
    xx=np.linspace(x[0],x[-1],2000)
    
    f=interp1d(x,y,kind='cubic')
    yint=copy(f(xx))
    p.plot(xx,yint,label='interp1d')
    
    
    f=InterpolatedUnivariateSpline(x,y)
    yIUS=f(xx)
    p.plot(xx,yIUS,label='InterpolatedUnivariateSpline')
    
    f=UnivariateSpline(x,y)
    yUS=f(xx)
    p.plot(xx,yUS,label='UnivariateSpline')
    
    spl, X=fit_spline(x,y)
    xspl,yspl=plot_spline(spl,X)
    p.plot(xspl,yspl,label='fit_spline')
    
    a,b,c,d=np.polyfit(x,y,3)
    ypl=a*(xx**3)+b*(xx**2)+c*xx+d
    p.plot(xx,ypl,label='polyfit')
    p.legend()
    
def test_polys(x,y, N=5):
    #p.figure()
    p.plot(x,y,'ro')
    xx=np.linspace(x[0],x[-1],10*len(x))
    for i in np.arange(N)+1:
        coeffs=np.polyfit(x,y,i)
        yy=np.zeros_like(xx)
        for j in range(len(coeffs)):
            yy+=coeffs[i-j]*(xx**j)
        p.plot(xx,yy,label='degree-%d' %i)
    #p.legend()
        
