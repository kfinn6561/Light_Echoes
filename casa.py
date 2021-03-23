'''
Created on Jul 28, 2015

@author: kieran
'''
import numpy as np
import pylab as p
from general_tools import pload,reduce_range,pdump,mat_x_vec
from scipy.interpolate import interp1d,UnivariateSpline
import spectral_tools as st
from copy import copy
import sys
import warnings
warnings.simplefilter('ignore', np.RankWarning)



def remove_continuum(x,y):
    N=13#want a 13 point spline
    lmin=2500.
    lmax=10000.
    logls=np.linspace(np.log(lmin),np.log(lmax),N+1)
    ls=np.exp(logls)#want the ls to be linearly spaced in log space
    xx=[]
    yy=[]
    for i in range(len(ls)-1):
        lo,hi=reduce_range(x,ls[i],ls[i+1],indices=True)
        if lo!=hi:
            xx.append(np.mean(x[lo:hi]))
            yy.append(np.mean(y[lo:hi]))
    
    
    cont=UnivariateSpline(xx,yy)(x)
    return(y-cont)
    
def scale(x,*args):#scale such that the average value is 1, er1 and er2 are errors that are scaled the same way
    offset=0.
    x=np.array(x)
    norm=np.mean(x)
    out=[(x-offset)/norm]
    for arg in args:
        out.append((arg-offset)/norm)
    return out

c=2.99792458e2#10^3 km/s
def lambdas_to_vs(lambdas,l0):
    R=(lambdas/l0)**2
    return c*(R-1)/(R+1)
def vs_to_lambdas(vs,l0):
    beta=vs/c
    return l0*np.sqrt((1+beta)/(1-beta))

def onoff(x):
    if x:
        return 'on'
    else:
        return 'off'

le_color='b'
casa_color='#CD0000'


onefigure=not ('multifig' in sys.argv) 
    
VERBOSE=('verbose' in sys.argv)
show_ranges=('ranges' in sys.argv)
full_range=('fullfit' in sys.argv)
smoothing=not ('unsmooth' in sys.argv)
raw=('raw' in sys.argv)
novel=('no_velocity'in sys.argv)
rem_cont=('remove_continuum' in sys.argv)
template_only=('template_only' in sys.argv)
show_plots=not ('no_plots' in sys.argv)


if 'minimum' in sys.argv:
    polynom='minimum'
elif 'cubic' in sys.argv:
    polynom='cubic_poly'
else:
    polynom='quadratic_poly'
if VERBOSE or novel:
    save_plotting=('plotting' in sys.argv)
else:
    save_plotting=not ('trial' in sys.argv)#use trial argument then data is not overwritten 
 

print 'mutifig is %s' %onoff(not onefigure)
print 'verbose is %s' %onoff(VERBOSE)
print 'show_ranges is %s' %onoff(show_ranges)
print 'full_range is %s' %onoff(full_range)
print 'smoothing is %s' %onoff(smoothing)
print 'save_plotting is %s' %onoff(save_plotting)
print 'raw is %s' %onoff(raw)
print 'no_velocity is %s' %onoff(novel)
print 'remove_continuum is %s' %onoff(rem_cont)
print 'template_only is %s' %onoff(template_only)
print 'using %s fitting' %polynom


st.VERBOSE=VERBOSE

#le_ranges={'He':(5545,5860),'H_alpha':(6310,6460),'Ca_III':(8200,8670),'Fe_II':(5000,5150)}

le_ranges={'He':(5545,5860),'H_alpha':(6200,6400),'Ca_III':(8200,8670),'Fe_II':(5000,5150)}


casa_ranges={'He':{'le2116':(5600,5850),'le3923':(5555,5880),'le2521':(5600,5850)},
             'H_alpha':{'le2116':(6230,6430),'le3923':(6285,6470),'le2521':(6230,6500)},
             'Ca_III':{'le2116':(8080,8625),'le3923':(8280,8425),'le2521':(8120,8550)},
             'Fe_II':{'le2116':(5090,5150),'le3923':(5025,5100),'le2521':(5060,5140)}}


'''
casa_ranges={'He':{'le2116':(5600,5850),'le3923':(5555,5880),'le2521':(5600,5850)},
             'H_alpha':{'le2116':(6230,6430),'le3923':(6285,6470),'le2521':(6125,6500)},
             'Ca_III':{'le2116':(8080,8625),'le3923':(8280,8425),'le2521':(8120,8550)},
             'Fe_II':{'le2116':(5090,5150),'le3923':(5025,5100),'le2521':(5060,5140)}}


casa_centres={'H_alpha':{'le2116':6320,'le2521':6300},'He':{'le2116':5730}}
le_centres={'H_alpha':{'le2116':6320,'le2521':6320,'le3923':6320},
            'He':{'le2116':5760,'le2521':5760,'le3923':5760}}
'''


casa2116=np.loadtxt("../spectra/casales/casa2116-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
casa2521=np.loadtxt("../spectra/casales/casa2521-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})
casa3923=np.loadtxt("../spectra/casales/casa3923-20090922.flm", dtype={'names': ('l', 'f', 'var'), 'formats': ('f4', 'f4', 'f4')})

raw_2116=np.loadtxt('../reducedspec/tyc2116-20090922.006-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
#raw_2521=np.loadtxt('../reducedspec/tyc2521-20090922.005-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})
raw_2521=pload('combined_2521.dat')
raw_3923=np.loadtxt('../reducedspec/tyc3923-20091023.503-br.flm', dtype={'names': ('l', 'f'), 'formats': ('f4', 'f4')})

raw_data={'le2116':raw_2116,'le2521':raw_2521,'le3923':raw_3923}

casa_data={'le2116':casa2116,'le2521':casa2521,'le3923':casa3923}
les=['le3923','le2521','le2116']

casa_functions={}
casa_sig_functions={}
raw_functions={}
for le in les:
    if rem_cont:
        raw_d=remove_continuum(raw_data[le]['l'],raw_data[le]['f'])
    else:
        raw_d=raw_data[le]['f']
    casa_functions[le]=interp1d(casa_data[le]['l'],casa_data[le]['f'])
    #casa_data[le]['var']/=2#####not sure if this is correct but otherwise it doesn't fit the data
    
    casa_sig_functions[le]=interp1d(casa_data[le]['l'],casa_data[le]['var'])
    #raw_functions[le]=interp1d(raw_data[le]['l'],raw_data[le]['f'])
    raw_functions[le]=interp1d(raw_data[le]['l'],raw_d)
    ''''p.figure()
    p.plot(raw_data[le]['l'],raw_rem_cont)
    p.title(le)'''


lines=pload('features_le.pkl')

line_names=['He','H_alpha','Ca_III']


'''
line_names=['H_alpha','Ca_III']
les=['le2521','le3923']
'''

line_plotnames={'He':'He 5876','H_alpha':r'$H\alpha$','Ca_III':'Ca II','Fe_II':'Fe 5169'}

cols=len(line_names)


if onefigure:
    f,ax=p.subplots(len(les),cols,sharex='col',sharey='row')
    #f,ax=p.subplots(len(les),cols,sharex='col',sharey='row')
    p.tight_layout(pad=3,h_pad=-1,w_pad=-1)
    for axx in ax.flatten():
        axx.minorticks_on()
    ax=ax.T
else:
    ax=[]
    figs=[]
    for i in range(cols):
        f,axx=p.subplots(len(les),1,sharex=True)
        ax.append(axx)
        figs.append(f)
    for axx in np.array(ax).flatten():
        axx.minorticks_on()

ax2=np.ones_like(ax)
for i in range(cols):
    for j in range(len(les)):
        ax2[i][j]=ax[i][j].twiny()

plot_out=[[0 for le in les] for i in range(cols)]
velo_out={line:{le:{} for le in les} for line in line_names}
for i in range(cols):
    line=lines[line_names[i]]
    print 'calculating for %s line' %line_names[i]
    for j in range(len(les)):
        print 'Dust sheet %d of %d (%s)' %(j+1,len(les),les[j])
        l0=line['l0']
        le_lambdas=line['lambdas']
        full_lambdas=casa_data[les[j]]['l']
        casa_lambdas=reduce_range(full_lambdas,le_lambdas[0],le_lambdas[-1])
        
        
        casa_mid=casa_functions[les[j]](casa_lambdas)
        casa_sig=casa_sig_functions[les[j]](casa_lambdas)
        casa_raw=raw_functions[les[j]](casa_lambdas)
        casa_min=casa_mid-casa_sig
        casa_max=casa_mid+casa_sig
        #casa_mid,casa_max,casa_min=remove_continuum(lambdas,casa_mid,casa_max,casa_min)
        casa_mid,casa_max,casa_min,casa_raw=scale(casa_mid,casa_max,casa_min,casa_raw)
        
        
        le_mid,le_min,le_max=line[les[j]]
        le_mid,le_min,le_max=scale(le_mid,le_min,le_max)
        le_sig=(le_max-le_min)/2.#get_min currently only takes a symmetric varience, not a min and max
        casa_sig=(casa_max-casa_min)/2.#the scaling process changes the varience. This gets the real ones
        
        
        le_X=le_lambdas
        casa_X=casa_lambdas
        
        
        if VERBOSE:
            p.figure()
            
        print 'Finding minimum for the template spectrum'
        try:
            le_range=le_ranges[line_names[i]][les[j]]
        except TypeError:
            le_range=le_ranges[line_names[i]]
        if not novel:
            le_x0,le_xmin,le_xmax=st.get_min(le_lambdas,le_mid,le_sig,fixed_range=le_range,
                                full_range=full_range,smoothing=False,fitting=polynom)
        else:
            le_x0=le_xmax=le_xmin=0
            
        if VERBOSE:
            p.plot(le_lambdas,le_mid,'k')
            p.plot(le_lambdas,le_mid,'ko')
            p.fill_between(le_lambdas,le_min,le_max,color='r',alpha=0.4)
            p.title('Template: line %s in %s' %(line_names[i],les[j]))
        
            p.figure()
        print 'finding minimum for the Cas A spectrum'
        casa_mean,casa_cov,SM=st.get_regression(lambdas_to_vs(casa_lambdas,l0),casa_raw,casa_sig)
        '''
        if line_names[i]=='He' and les[j]=='le2521':
            kdjfghdskhg
        '''
        if novel or template_only:
            casa_x0=casa_xmin=casa_xmax=0
        else:
            casa_x0,casa_xmin,casa_xmax=st.get_min(casa_lambdas,casa_mean,casa_cov,fixed_range=casa_ranges[line_names[i]][les[j]],
                                full_range=False,smoothing=smoothing,fitting=polynom,smoothing_matrix=SM)#change fitting=polynom, full_range=full_range
    
        
        '''Rebinning'''
        casa_X,casa_mean,casa_cov,A=st.rebin(casa_X,le_X,casa_mean,casa_cov,return_A=True)
        SM=A*SM*(A.I)#the smoothing matrix for the rebinned spectrum
        casa_lambdas=casa_X#todo, see why both casa_X and casa_lambdas exist and get rid of the unecessary one
        
        
        
        
        casa_newvar=np.diagonal(casa_cov)
        casa_newsig=np.sqrt(casa_newvar)
        #casa_mean=casa_mean+1.#offset
        casa_min=casa_mean-casa_newsig
        casa_max=casa_mean+casa_newsig    
            
        if VERBOSE:
            p.plot(casa_lambdas,casa_mean,'k')
            p.plot(casa_lambdas,casa_mean,'ko')
            #p.plot(casa_lambdas,casa_mid,'b',lw=2)
            p.fill_between(casa_lambdas,casa_min,casa_max,color='r',alpha=0.4)
            p.title('Cas A: line %s in %s' %(line_names[i],les[j]))
        
        velocities=np.append(lambdas_to_vs(le_lambdas, l0),lambdas_to_vs(le_lambdas, l0))
        x=min(velocities)
        minv=copy(x-x%5+5)
        velo_labels=np.arange(minv,max(velocities),5)
        velo_positions=vs_to_lambdas(velo_labels, l0)
        ax2[i][j].set_xticks(velo_positions)
        if j==0:
            ax2[i][j].set_xticklabels([str(int(np.round(vl))) for vl in velo_labels])
            if i==1:
                ax2[i][j].set_xlabel('Velocity (10^3 km/s)')
        else:
            ax2[i][j].set_xticklabels([])
        
        ax[i][j].plot(le_X,le_mid,color=le_color,label='Type IIb template')
        ax[i][j].fill_between(le_X,le_min,le_max,color=le_color,alpha=0.4)
        ax[i][j].axvline(le_x0,color=le_color)
        ax[i][j].axvspan(le_xmin,le_xmax,facecolor=le_color,alpha=0.4)
        
        ax[i][j].plot(casa_X,casa_mean,color=casa_color,label='Cas A spectrum')
        ax[i][j].fill_between(casa_X,casa_min,casa_max,color=casa_color,alpha=0.3)
        ax[i][j].axvline(casa_x0,color=casa_color)
        ax[i][j].axvspan(casa_xmin,casa_xmax,facecolor=casa_color,alpha=0.3)
        plot_out[i][j]=[[velo_labels,velo_positions],
                        [le_X,le_mid,le_min,le_max,le_x0,le_xmin,le_xmax],
                        [casa_X,casa_mean,casa_min,casa_max,casa_x0,casa_xmin,casa_xmax]]
        
        casa_xv,casa_xvmin,casa_xvmax=lambdas_to_vs(np.array([casa_x0,casa_xmin,casa_xmax]), l0)
        velo_out[line_names[i]][les[j]]['Cas A']=[casa_xv,casa_xvmax,casa_xvmin]
        
        le_xv,le_xvmin,le_xvmax=lambdas_to_vs(np.array([le_x0,le_xmin,le_xmax]), l0)
        velo_out[line_names[i]][les[j]]['Template']=[le_xv,le_xvmax,le_xvmin]
        
        if VERBOSE or show_ranges:
            casa_left,casa_right=casa_ranges[line_names[i]][les[j]]
            le_left,le_right=le_range
            ax[i][j].axvline(casa_left,color=casa_color,ls='--')
            ax[i][j].axvline(casa_right,color=casa_color,ls='--')
            ax[i][j].axvline(le_left,color=le_color,ls='--')
            ax[i][j].axvline(le_right,color=le_color,ls='--')
       
        ax[i][j].set_xlim(le_X[0],le_X[-1])
        ax[i][j].set_ylim(0,2)
        if j!=0:
            ax[i][j].set_yticks([0,0.5,1,1.5])
        ax[i][j].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax2[i][j].set_xlim(ax[i][j].get_xlim())
        if not onefigure:
            ax[i][j].set_ylabel(les[j])
       
    if onefigure:
        ax[i][0].text(.05,0.05,line_plotnames[line_names[i]],horizontalalignment='left',transform=ax[i][0].transAxes)
    else:
        ax[i][0].set_title(line_names[i])
        ax[i][0].legend(loc='upper left')
        figs[i].savefig('casa_%s.png' %line_names[i])

if onefigure:    
    for j in range(len(les)):
        ax[0][j].set_ylabel(les[j])
    ax[-1][0].legend(loc='upper left')
    ax[1][-1].set_xlabel('Wavelength (A)')
    p.savefig('casa.png')
    
if save_plotting and (not VERBOSE):
    pdump([les,line_names,plot_out],'plot_data/casa.pkl')
    pdump(velo_out,'plot_data/casa_velocity_list.pkl')

if show_plots:
    p.show()       