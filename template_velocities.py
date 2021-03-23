'''
Created on Jul 28, 2015

@author: kieran
'''
import numpy as np
import pylab as p
from general_tools import pload,reduce_range
from scipy.interpolate import interp1d
import spectral_tools as st
from copy import copy
import sys
import warnings
warnings.simplefilter('ignore', np.RankWarning)

onefigure=True 
velocity=not ('wavelength' in sys.argv)
show_ranges=('ranges' in sys.argv)
full_range=False
st.RECALC=False
smoothing= True
if 'quadratic' in sys.argv:
    polynom='quadratic_poly'
else:
    polynom='cubic_poly'

 

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
    
def plot_min(lambdas,flux,var,range,l0,ax,label,color='b',style='-'):
    x0,err=st.get_min(lambdas,flux,var,fixed_range=range,
                                full_range=full_range,smoothing=smoothing,fitting=polynom)
    if velocity:
        x0,err=lambdas_to_vs(np.array([x0,err]),l0)
    ax.axvline(x0,color=color,ls=style,label=label)
    

le_ranges={'He':(5545,5860),'H_alpha':(6100,6490),'Ca_III':(8200,8670)}
le_color='r'
les=np.sort(['le2116','le2521','le3923'])

lines=pload('extreme.pkl')

line_names=['He','H_alpha','Ca_III']

cols=len(line_names)
f,ax=p.subplots(len(les),cols,sharex='col')
ax=ax.T
color_cycle = ax[0][0]._get_lines.color_cycle


for i in range(cols):
    line=lines[line_names[i]]
    print 'calculating for %s line' %line_names[i]
    for j in range(len(les)):
        print 'Dust sheet %d of %d' %(j+1,len(les))
        l0=line['l0']
        le_lambdas=line['lambdas']        
        
        le_mid,le_min,le_max,extreme_min,extreme_max=line[les[j]]
        le_mid,le_min,le_max,extreme_min,extreme_max=scale(le_mid,le_min,le_max,extreme_min,extreme_max)
        le_var=(le_max-le_min)/2.#get_min currently only takes a symmetric varience, not a min and max
        
        if velocity:
            le_X=lambdas_to_vs(le_lambdas,l0)
        else:
            le_X=le_lambdas
        
        ax[i][j].plot(le_X,le_mid,color=le_color)
        ax[i][j].fill_between(le_X,le_min,le_max,color=le_color,alpha=0.6)
        ax[i][j].fill_between(le_X,extreme_min,extreme_max,color=le_color,alpha=0.3)
        
        print 'Finding minimum for the mean spectrum'   
        plot_min(le_lambdas,le_mid,le_var,le_ranges[line_names[i]],l0,ax[i][j],'mean',color='k')
        plot_min(le_lambdas,le_min,le_var,le_ranges[line_names[i]],l0,ax[i][j],'-sd',color='r')
        plot_min(le_lambdas,le_max,le_var,le_ranges[line_names[i]],l0,ax[i][j],'+sd',color='m')
        plot_min(le_lambdas,extreme_max,le_var,le_ranges[line_names[i]],l0,ax[i][j],'max',color='b')
        plot_min(le_lambdas,extreme_min,le_var,le_ranges[line_names[i]],l0,ax[i][j],'min',color='g')
    
        ax[i][j].set_ylim(0,2)
        if not velocity:
            ax[i][j].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        
        if not onefigure:
            ax[i][j].set_ylabel(les[j])
       
    if onefigure:    
        ax[i][0].set_title(line_names[i])
        #ax[i][1].set_ylim(0.5,1.5)
    else:
        ax[i][0].set_title(line_names[i])
        ax[i][0].legend(loc='upper left')
    if velocity:
        ax[i][-1].set_xlabel('Velocity (10^3 km/s)')
    else:
        ax[i][-1].set_xlabel('Wavelength (A)')

if onefigure:    
    for j in range(len(les)):
        ax[0][j].set_ylabel(les[j])
    ax[1][0].legend(loc='upper center',ncol=5,fontsize=12)

p.show()       
        
        
        
        