import pylab as p
import numpy as np
from general_tools import pload,reduce_range,pdump
from scipy.interpolate import interp1d
from copy import copy

import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})#increase font size
mpl.rcParams['figure.max_open_warning']=False
le_plotnames={'le3923':'NE','le2521':'NW','le2116':'SW','unmodified':'w(t)=1'}

c=2.99792458e2#10^3 km/s
Ns=[0.01,0.1,0.4,1.,10.]
ls=[0.1,1.,2.,3.,4.]

def get_regression(x,y,var,l,N):
    '''uses regression to estimate the mean and the variance of the underlying distribution
    at the same positions as measured. Assumes gaussian kernel of width l'''
    
    #vexp=l_mult*0.00333
    #l=vexp*np.mean(x)#width of the kernel
    k=lambda a,b: N*np.exp(-0.5*((a-b)/l)**2)#gaussian kernel
    
    x=np.array(x)
    y=np.array(y)
    offset=copy(np.mean(y))
    y=y-offset
    y=np.matrix(y)
    var=np.diag(var)
    I=np.ones_like(x)
    Xi=np.outer(I,x)
    Xj=np.outer(x,I)
    
    K=k(Xi,Xj)#nxn matrix of covarience values
    K_er=np.matrix(copy(K)+var).I#this is the bottleneck, inverting the matrix
    K_xx=np.matrix(K)
    mean=K_xx*K_er*y.T
    cov=K_xx-K_xx*K_er*K_xx+var
    return (np.array(mean).T[0]+offset,np.diagonal(np.array(cov)))

def scale(x,*args):#scale such that the average value is 1, er1 and er2 are errors that are scaled the same way
    x=np.array(x)
    offset=0.
    norm=np.mean(x)
    out=[(x-offset)/norm]
    for arg in args:
        out.append((arg-offset)/norm)
    return out

def lambdas_to_vs(lambdas,l0):
    R=(lambdas/l0)**2
    return c*(R-1)/(R+1)

casa2116=np.loadtxt("../spectra/casales/casa2116-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
casa2521=np.loadtxt("../spectra/casales/casa2521-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
casa3923=np.loadtxt("../spectra/casales/casa3923-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})

raw_2116=np.loadtxt('../reducedspec/tyc2116-20090922.006-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
raw_2521=np.loadtxt('../reducedspec/tyc2521-20090922.005-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
raw_3923=np.loadtxt('../reducedspec/tyc3923-20091023.503-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})

raw_data={'le2116':raw_2116,'le2521':raw_2521,'le3923':raw_3923}

casa_data={'le2116':casa2116,'le2521':casa2521,'le3923':casa3923}
les=['le3923','le2521','le2116']


try:
    data=pload('regression_params.pkl')
except:
    data={}

casa_functions={}
casa_var_functions={}
raw_functions={}
for le in les:
    casa_functions[le]=interp1d(casa_data[le]['l'],casa_data[le]['f'])
    casa_var_functions[le]=interp1d(casa_data[le]['l'],casa_data[le]['var'])
    raw_functions[le]=interp1d(raw_data[le]['l'],raw_data[le]['f'])

lines=pload('features_le.pkl')

line_names=['He','H_alpha','Ca_III']
line_plotnames={'He':'He 5876','H_alpha':r'$H\alpha$','Ca_III':'Ca II','Fe_II':'Fe 5169'}

for ln in line_names:
    line=lines[ln]
    print 'calculating for %s line' %ln
    for le in les:
        print le
        l0=line['l0']
        le_lambdas=line['lambdas']
        full_lambdas=casa_data[le]['l']
        casa_lambdas=reduce_range(full_lambdas,le_lambdas[0],le_lambdas[-1])
        casa_mid=casa_functions[le](casa_lambdas)
        casa_var=casa_var_functions[le](casa_lambdas)
        casa_min=casa_mid-casa_var
        casa_max=casa_mid+casa_var
        casa_raw=raw_functions[le](casa_lambdas)
        casa_mid,casa_max,casa_min,casa_raw=scale(casa_mid,casa_max,casa_min,casa_raw)
        casa_var=(casa_max-casa_min)/2.#the scaling process changes the varience. This gets the real ones
                

        
       
        f,ax=p.subplots(len(ls),len(Ns),sharex='col',sharey='row',figsize=(18,12))
        p.tight_layout(pad=3,h_pad=-1,w_pad=-1)
        ax=ax.T
        for i in range(len(Ns)):
            ax[i][0].text(.05,0.82,'N=%g' %Ns[i],horizontalalignment='left',transform=ax[i][0].transAxes, fontsize=24)
            ax[i][-1].set_xlabel(r'Wavelength ($10^3\AA$)') 
            for j in range(len(ls)):
                print '%d of %d' %(1+i*len(ls)+j,len(Ns)*len(ls))
                label='regression_%s_%s_N=%g_l=%g' %(ln,le,Ns[i],ls[j])
                try:
                    mean,cov=data[label]
                except KeyError:
                    mean,cov=get_regression(lambdas_to_vs(casa_lambdas,l0),casa_raw,casa_var,ls[j],Ns[i])
                    data[label]=[mean,cov]
                ax[i][j].plot(casa_lambdas/1000,mean,color='k')
                ax[i][j].fill_between(casa_lambdas/1000,mean-cov,mean+cov,color='k',alpha=0.4)
                ax[i][j].errorbar(casa_lambdas/1000,casa_mid,color='b',yerr=casa_var,alpha=0.2)
                
                ax[i][j].set_xlim(le_lambdas[0]/1000,le_lambdas[-1]/1000)
                ax[i][j].set_ylim(0,1.99)
                
         
        for j in range(len(ls)):           
            ax[0][j].set_ylabel(r'Scaled $f_\lambda$')
            ax[0][j].text(.05,0.1,'l=%g'%ls[j],horizontalalignment='left',transform=ax[0][j].transAxes, fontsize=24)
        f.suptitle('%s line in light echo %s' %(line_plotnames[ln],le_plotnames[le]))
        f.savefig('regression_plots/regression_%s_%s.png' %(ln,le))
pdump(data,'regression_params.pkl')
p.show()
print('done')