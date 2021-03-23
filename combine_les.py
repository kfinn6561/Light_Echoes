'''
Created on Jul 28, 2015

@author: kieran
'''
import numpy as np
import pylab as p
from general_tools import pload,reduce_range,pdump
from scipy.interpolate import interp1d
from spectral_tools import smooth_spectrum, rebin
from copy import copy

def scale(x,*args):#scale such that the average value is 1, er1 and er2 are errors that are scaled the same way
    offset=0.
    x=np.array(x)
    norm=np.mean(x)
    out=[(x-offset)/norm]
    for arg in args:
        out.append((arg-offset)/norm)
    return out

gr=(np.sqrt(5)-1)/2
def gss(f,a,b,tol=1e-5):
    '''
    golden section search
    to find the minimum of f on [a,b]
    f: a strictly unimodal function on [a,b]

    example:
    >>> f=lambda x:(x-2)**2
    >>> x=gss(f,1,5)
    >>> x
    2.000009644875678

    '''
    c=b-gr*(b-a)
    d=a+gr*(b-a)
    while abs(c-d)>tol:       
        fc=f(c);fd=f(d)
        if fc<fd:
            b=d
            d=c  #fd=fc;fc=f(c)
            c=b-gr*(b-a)
        else:
            a=c
            c=d  #fc=fd;fd=f(d)
            d=a+gr*(b-a)
    return (b+a)/2

raw_1=np.loadtxt('../reducedspec/tyc2521-20090922.005-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
raw_2=np.loadtxt('../reducedspec/tyc2521mmt-20090921.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})



smooth_1=np.loadtxt("../spectra/casales/casa2521-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
smooth_2=np.loadtxt("../spectra/casales/casa2521-20090921.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})


smooth_2['var']*=1e-9#I'm not sure if this is the correct factor but something is wrong about the errors in this file.
#I should probably not use these errors for anything important

reduced_lambdas=reduce_range(smooth_2['l'],3500,8350)
smooth_function_1=interp1d(smooth_1['l'],smooth_1['f'])
compare_1=smooth_function_1(reduced_lambdas)

smooth_function_2=interp1d(smooth_2['l'],smooth_2['f'])
raw_function_2=interp1d(raw_2['l'],raw_2['f'])
var_function_2=interp1d(smooth_2['l'],smooth_2['var'])
compare_2=smooth_function_2(reduced_lambdas)

p.plot(raw_1['l'],raw_1['f'])
p.plot(smooth_1['l'],smooth_1['f'])
p.fill_between(smooth_1['l'],smooth_1['f']-smooth_1['var'],smooth_1['f']+smooth_1['var'],alpha=0.5)

p.figure()
p.plot(raw_2['l'],raw_2['f'])
p.plot(smooth_2['l'],smooth_2['f'])
p.fill_between(smooth_2['l'],smooth_2['f']-smooth_2['var'],smooth_2['f']+smooth_2['var'],alpha=0.5)

def compare_spectra(k):#k is the multiplicative constant for spectrum_2
    f2=compare_2*k
    out=(f2-compare_1)**2
    return sum(out)

p.figure()
#x=np.logspace(14,18,1000)
x=np.linspace(1,3,100)*1e15
y=[compare_spectra(k) for k in x]
p.plot(x,y)
#p.loglog()

optimal_k=gss(compare_spectra,1e15,3e15,1e10)
p.plot(optimal_k,compare_spectra(optimal_k),'ko')


p.figure()
p.plot(smooth_1['l'],smooth_1['f'])
p.plot(smooth_2['l'],smooth_2['f']*optimal_k)

out_f=[]
out_var=[]
lambdas=copy(raw_1['l'][1::2])

cov=np.matrix(np.diag(smooth_1['var']))
raw1_rb={}
raw1_rb['l'],raw1_rb['f'],var=rebin(raw_1['l'],lambdas,raw_1['f'],cov)
raw1_rb['var']=np.diagonal(var)#may need to use full covariance. Shouldn't be necessary though


for i in range(len(raw1_rb['l'])):
    l=raw1_rb['l'][i]
    if l in raw_2['l']:
        out_f.append(0.5*(raw1_rb['f'][i]+optimal_k*raw_function_2(l)))
        out_var.append(0.5*np.sqrt(raw1_rb['var'][i]**2+(var_function_2(l)*optimal_k)**2))
    else:
        out_f.append(raw1_rb['f'][i])
        out_var.append(raw1_rb['var'][i])
        
lambdas=np.array(lambdas)
out_f=np.array(out_f)
out_var=np.array(out_var)
out_smooth=smooth_spectrum(lambdas,out_f,out_var)        
        
p.figure()
p.plot(lambdas,out_f)
p.fill_between(lambdas,out_smooth+out_var,out_smooth-out_var,alpha=0.4)

pdump({'l':lambdas,'f':out_f,'var':out_var},'combined_2521.dat')

p.show()