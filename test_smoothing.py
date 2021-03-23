'''
Created on Aug 26, 2015

@author: kieran
'''
import pylab as p
import numpy as np
from spectral_tools import smooth_spectrum
from general_tools import reduce_range,progress_bar

Nscatters=100

real=np.loadtxt('../reducedspec/tyc2116-20090922.006-smooth-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})

lo,hi=reduce_range(real['l'],7000.,8000., indices=True)
lambdas=real['l'][lo:hi]
mean=real['f'][lo:hi]
mean=mean/max(mean)
mean+=np.abs(np.min(mean))

var=np.sqrt(np.abs(mean))#poisson noise

p.figure()
p.plot(lambdas,mean)
p.fill_between(lambdas,mean-var,mean+var,alpha=0.4)

raw=np.random.normal(mean,var)
new_var=np.sqrt(np.abs(raw))

p.plot(lambdas, raw,'ko')

p.figure()

p.plot(lambdas,mean)
p.fill_between(lambdas,mean-var,mean+var,alpha=0.4)
for i in range(Nscatters):
    progress_bar(i,Nscatters)
    tflux=np.random.normal(raw,new_var)
    y=smooth_spectrum(lambdas, tflux, new_var)
    p.plot(lambdas,y,'rx')
#p.title('raw')

print '\n'
#p.figure()
smooth=smooth_spectrum(lambdas,raw,var)
#p.plot(lambdas,mean)
#p.fill_between(lambdas,mean-var,mean+var,alpha=0.4)
for i in range(Nscatters):
    progress_bar(i,Nscatters)
    tflux=np.random.normal(smooth,new_var)
    y=smooth_spectrum(lambdas, tflux, new_var)
    p.plot(lambdas,y,'bx')
#p.title('smooth')

p.figure()

p.plot(lambdas,mean,color='b')
p.fill_between(lambdas,mean-var,mean+var,alpha=0.4,color='b')
p.plot(lambdas,raw,color='r')
p.fill_between(lambdas,raw-var,raw+var,alpha=0.4,color='r')
p.title('raw')

p.figure()
p.plot(lambdas,mean,color='b')
p.fill_between(lambdas,mean-var,mean+var,alpha=0.4,color='b')
p.plot(lambdas,smooth,color='r')
p.fill_between(lambdas,smooth-var,smooth+var,alpha=0.4,color='r')
p.title('smooth')
p.show()


