'''
Created on Aug 21, 2015

@author: kieran
'''
import pylab as p
import numpy as np
from general_tools import reduce_range,increase_sampling,pdump

le='le2116'
line='H_alpha'

def scale(x,er1,er2):#scale such that the average value is 1, er1 and er2 are errors that are scaled the same way
    x=np.array(x)
    offset=0.
    norm=np.mean(x)
    return ((x-offset)/norm,(er1-offset)/norm,(er2-offset)/norm)


ranges={'He':{'le2116':(5600,5850),'le3923':(5555,5880),'le2521':(5600,5850)},
             'H_alpha':{'le2116':(6230,6430),'le3923':(6285,6470),'le2521':(6125,6500)},
             'Ca_III':{'le2116':(8080,8625),'le3923':(8280,8425),'le2521':(8120,8550)}}


casa2116=np.loadtxt("../spectra/casales/casa2116-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
casa2521=np.loadtxt("../spectra/casales/casa2521-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
casa3923=np.loadtxt("../spectra/casales/casa3923-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})


data={'le2116':casa2116,'le2521':casa2521,'le3923':casa3923}

lambdas=data[le]['l']
mid=data[le]['f']
var=data[le]['var']

left,right=ranges[line][le]
lo,hi=reduce_range(lambdas,left,right,indices=True)
lambdas=lambdas[lo:hi]
mid=mid[lo:hi]
var=var[lo:hi]
minn=mid-var
maxx=mid+var
mid,minn,maxx=scale(mid,minn,maxx)
var=(maxx-minn)/2.


tflux=np.random.normal(mid,var)

p.plot(lambdas,mid,color='k',lw=2)
p.fill_between(lambdas,minn,maxx,alpha=0.4)
p.plot(lambdas,tflux,'rx')


a,b,c=np.polyfit(lambdas,tflux,2)
print a,b,c
x=increase_sampling(lambdas,10)
y=a*(x**2)+b*x+c
p.plot(x,y)

pdump([lambdas,tflux],'random_sample.dat')

p.show()




 