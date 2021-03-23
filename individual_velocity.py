'''
Created on Sep 16, 2015

@author: kieran
'''
import numpy as np
import pylab as p
from general_tools import pload,reduce_range,pdump
from scipy.interpolate import interp1d
import spectral_tools as st
from copy import copy
import sys
import warnings
from colour_tools import kelly_colors
warnings.simplefilter('ignore', np.RankWarning)

import matplotlib as mpl
mpl.rcParams.update({'font.size': 18})#increase font size
mpl.rcParams['xtick.major.size'] = 20
mpl.rcParams['xtick.major.width'] = 4
mpl.rcParams['xtick.minor.size'] = 10
mpl.rcParams['xtick.minor.width'] = 2#tick size
mpl.rcParams['ytick.major.size'] = 20
mpl.rcParams['ytick.major.width'] = 4
mpl.rcParams['ytick.minor.size'] = 10
mpl.rcParams['ytick.minor.width'] = 2#tick size


def remove_continuum(x,y,er1,er2):
    f=interp1d(x,y)
    xx=np.linspace(x[0],x[-1],13)#want a 13 point spline
    cont=interp1d(xx,f(xx),kind='cubic')(x)
    return(y/cont,er1/cont,er2/cont)
    
def scale(x,*args):#scale such that the average value is 1, er1 and er2 are errors that are scaled the same way
    x=np.array(x)
    offset=0.
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
smoothing=('smooth' in sys.argv)
raw=('raw' in sys.argv)
if 'minimum' in sys.argv:
    polynom='minimum'
elif 'cubic' in sys.argv:
    polynom='cubic_poly'
else:
    polynom='quadratic_poly'
save_plotting=not ('trial' in sys.argv)#use trial argument then data is not overwritten   
show_plots=not ('no_plots' in sys.argv)
 

print 'verbose is %s' %onoff(VERBOSE)
print 'show_ranges is %s' %onoff(show_ranges)
print 'full_range is %s' %onoff(full_range)
print 'smoothing is %s' %onoff(smoothing)
print 'save_plotting is %s' %onoff(save_plotting)
print 'raw is %s' %onoff(raw)
print 'using %s fitting' %polynom

st.VERBOSE=VERBOSE

le_ranges={'He':(5545,5860),'H_alpha':(6260,6410),'Ca_III':(8200,8670),'Fe_II':(5000,5150)}

casa_ranges={'He':{'le2116':(5600,5850),'le3923':(5555,5880),'le2521':(5600,5850)},
             'H_alpha':{'le2116':(6230,6430),'le3923':(6285,6470),'le2521':(6125,6500)},
             'Ca_III':{'le2116':(8080,8625),'le3923':(8280,8425),'le2521':(8120,8550)},
             'Fe_II':{'le2116':(5090,5150),'le3923':(5120,5200),'le2521':(5060,5140)}}

'''
casa_centres={'H_alpha':{'le2116':6320,'le2521':6300},'He':{'le2116':5730}}
le_centres={'H_alpha':{'le2116':6320,'le2521':6320,'le3923':6320},
            'He':{'le2116':5760,'le2521':5760,'le3923':5760}}
'''
les=['le3923','le2521','le2116']

sn_colours={}
colours=iter(kelly_colors)



lines=pload('individual.pkl')

line_names=['He','H_alpha','Ca_III']
line_plotnames={'He':'He 5876','H_alpha':r'$H\alpha$','Ca_III':'Ca II','Fe_II':'Fe 5169'}


cols=len(line_names)


if onefigure:
    #f,ax=p.subplots(len(les),cols,sharex='col',sharey='row',figsize=(17,25))
    f,ax=p.subplots(len(les),cols,sharex='col',sharey='row')
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
        print 'Dust sheet %d of %d' %(j+1,len(les))
        l0=line['l0']
        temp_lambdas=line[line.keys()[0]]['lambdas']
        #full_lambdas=casa_data[les[j]]['l']
        
        plot_sne=[]
        if line_names[i]=='Ca_III':
            temp_sne=['sn2003bg','sn2011ei']
        else:
            temp_sne=line.keys()
        for sn in temp_sne:
            if sn=='l0':
                continue
            try:
                le_color=sn_colours[sn]
            except KeyError:
                sn_colours[sn]=next(colours)
                le_color=sn_colours[sn]
            
            le_mid,le_min,le_max=line[sn][les[j]]
            le_mid,le_min,le_max=scale(le_mid,le_min,le_max)
            le_sig=(le_max-le_min)/2.#get_min currently only takes a symmetric varience, not a min and max
            
            
            le_lambdas=line[sn]['lambdas']
            le_X=le_lambdas
        
        
            if VERBOSE:
                p.figure()
                
            print 'Finding minimum for %s' %sn   
            le_x0,le_xmin,le_xmax=st.get_min(le_lambdas,le_mid,le_sig,fixed_range=le_ranges[line_names[i]],
                                    full_range=full_range,smoothing=smoothing,fitting=polynom)
                
            if VERBOSE:
                p.plot(le_lambdas,le_mid,'k')
                p.plot(le_lambdas,le_mid,'ko')
                p.fill_between(le_lambdas,le_min,le_max,color='r',alpha=0.4)
                p.title('%s: line %s in %s' %(sn,line_names[i],les[j]))
            
                p.figure()
            ax[i][j].plot(le_X,le_mid,color=le_color,label=sn)
            ax[i][j].fill_between(le_X,le_min,le_max,color=le_color,alpha=0.4)
            ax[i][j].axvline(le_x0,color=le_color)
            ax[i][j].axvspan(le_xmin,le_xmax,facecolor=le_color,alpha=0.4)
            plot_sne.append([sn,le_X,le_mid,le_min,le_max,le_x0,le_xmin,le_xmax])
            
            
            le_xv,le_xvmin,le_xvmax=lambdas_to_vs(np.array([le_x0,le_xmin,le_xmax]), l0)
            velo_out[line_names[i]][les[j]][sn]=[le_xv,le_xvmax,le_xvmin]
          
        
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

        
        if VERBOSE or show_ranges:
            le_left,le_right=le_ranges[line_names[i]]
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
        
        plot_out[i][j]=[[velo_labels,velo_positions],plot_sne]
       
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
    pdump([les,line_names,plot_out],'plot_data/individual_velocity.pkl')
    pdump(velo_out,'plot_data/individual_velocity_list.pkl')

if show_plots:
    p.show()    
        
        
        
        