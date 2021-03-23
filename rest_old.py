'''
Created on June 2nd, 2015

RUN SOURCE EXPORT.SOURCEME BEFORE STARTING


@author: kieran
'''
import pylab as p
import numpy as np
from functions import *
from echo_tools_rect import *
from general_tools import progress_bar,pload,pdump
from copy import copy

'''use the following units
time: years
length: light years
angles: arcsecs (for the most part)
'''

sl=-0.1*365
sr=0.3*365
T=0.#time of measurement. This is arbitrary, but fixes the origin of time.
source=Source(triangle(sl,sr),sl,sr,T-300.)#paper says 300y old event 10000ly away but then wouldn't have seen it yet?
Nres=50#resolution of the window functions etc

def calc_profiles(dust_sheets,plotname):
    print'\n\n'
    '''LE Profile'''
    print 'calculating profiles'
    rho_arcs=np.linspace(-7,5,Nres)
    fname='window_data/'+plotname+'_profile.dat'
    try:
        profile=pload(fname)
        print 'profile loaded from file'
    except IOError:
        profile=[[] for dust in dust_sheets]
        '''cannot do as a vector due to the way integrate has been programmed.
        If this proves too slow might want to look into it'''
        for i in range(len(dust_sheets)):
            print '\n'
            print 'Dust sheet %d of %d' %(i+1,len(dust_sheets))
            dust_sheets[i].fill_flux_array(source)
            profile_function=dust_sheets[i].profile()
            profile[i]=profile_function(rho_arcs)
        pdump(profile,fname)
    return [profile,rho_arcs]

def calc_windows(dust_sheets,plotname):
    '''Window function'''
    print '\n'
    print 'calculating the window functions'
    window_left=-150.#days
    window_right=200.
    slit=Slit(0.,1.)
    ts=np.linspace(window_left,window_right,Nres)
    fname='window_data/'+plotname+'_windows.dat'
    try:
        windows=pload(fname)
        print 'windows data loaded from file'
    except IOError:
        windows=[[] for dust in dust_sheets]
        for i in range(len(dust_sheets)):
            print '\n'
            print 'Dust sheet %d of %d' %(i+1,len(dust_sheets))
            window_function=dust_sheets[i].window(slit)
            windows[i]=normalise(window_function(ts))
        pdump(windows,fname)
    return [windows,ts]

def plot_dust(dust_sheets,variables,plotname,rows):    
    profile,rho_arcs=calc_profiles(dust_sheets,plotname)
    windows,ts=calc_windows(dust_sheets,plotname)
    window_left=min(ts)
    window_right=max(ts)
    
    '''convolution'''
    print '\n'
    print 'convolving'    
    light_curves=[]
    scs=normalise(source(ts))
    for i in range(len(dust_sheets)):
            light_curves.append(normalise(windows[i])*scs)
    
    ns=range(len(dust_sheets))
    N=int(len(dust_sheets)/rows)
    row_ns=[ns[i*N:(i+1)*N] for i in range(rows-1)]
    row_ns.append(ns[(rows-1)*N:])#there may be some left over after splitting them up into rows
            
    f,ax=p.subplots(rows,3,sharex='col',sharey='row')
    for row in range(rows):
        for i in row_ns[row]:
            ax[row][0].plot(rho_arcs,normalise(profile[i]),label=variables[i])
            ax[row][1].plot(ts,normalise(windows[i]),label=variables[i])
            ax[row][2].plot(ts,light_curves[i],label=variables[i])
        ax[row][2].plot(ts, normalise(source(ts)),'k--',label='light curve')
        ax[row][2].legend()
        ax[row][0].set_ylim(-0.1,1.1)
    ax[0][1].set_xlim(window_left,window_right)
    ax[0][2].set_xlim(-100,150)
    
    
    fig=p.gcf()
    fig.set_size_inches(20,20)
    fig.suptitle(plotname)
    p.savefig(plotname)




VARIABLES=[]
DUST_SHEETS=[]
PLOTNAMES=[]
ROWS=[]


dust_sheets=[]
variables=[0.,1./16,2./16,3./16,4./16]
rows=5
for centre in variables:
    for sig in [0.1,0.2]:
        dust_sheets.append(dust_sheet(sig,300.,0,45.,10000.,profile=lambda x:get_tophat(x,x*centre)))
variables=variables*2
plotname='kieran_off_centre'
DUST_SHEETS.append(copy(dust_sheets))
VARIABLES.append(copy(variables))
PLOTNAMES.append(copy(plotname))
ROWS.append(copy(rows))



dust_sheets=[]  
variables=[0.001,0.01,0.1,0.2,0.3,0.5]
rows=3
for sig in variables:
    dust_sheets.append(dust_sheet(sig,300.,0,45.,10000.))
plotname='kieran_widths'
DUST_SHEETS.append(copy(dust_sheets))
VARIABLES.append(copy(variables))
PLOTNAMES.append(copy(plotname))
ROWS.append(copy(rows))


dust_sheets=[]
variables=[70.,45.,0.,-30.]
rows=2
for alpha in variables:
    dust_sheets.append(dust_sheet(0.1,300.,0,alpha,10000.))
plotname='kieran_angles'
DUST_SHEETS.append(copy(dust_sheets))
VARIABLES.append(copy(variables))
PLOTNAMES.append(copy(plotname))
ROWS.append(copy(rows))


dust_sheets=[]
variables=[0,0.4,1.,2.]
sigmas=[0.001,0.1,0.3]
for sig in sigmas:
    for psf_sigma in variables:
        dust_sheets.append(dust_sheet(sig,300.,0,45.,10000.,psf=psf_sigma))
variables=variables*3
rows=3
plotname='kieran_seeing'
DUST_SHEETS.append(copy(dust_sheets))
VARIABLES.append(copy(variables))
PLOTNAMES.append(copy(plotname))
ROWS.append(copy(rows))


dust_sheets=[]
variables=[0.,-1.1,-2.2]
sigmas=[0.001,0.1,0.3]
for sig in sigmas:
    for ofst in variables:
        dust_sheets.append(dust_sheet(sig,300.,0,45.,10000.,offset=ofst))
variables=variables*3
rows=3
plotname='kieran_offsets'
DUST_SHEETS.append(copy(dust_sheets))
VARIABLES.append(copy(variables))
PLOTNAMES.append(copy(plotname))
ROWS.append(copy(rows))


dust_sheets=[]
variables=[0.,-1.1,-2.2]
sigmas=[0.001,0.1,0.3]
for sig in sigmas:
    for ofst in variables:
        dust_sheets.append(dust_sheet(sig,300.,0,45.,10000.,psf=1.,offset=ofst))
variables=variables*3
rows=3
plotname='kieran_offset_seeing'
DUST_SHEETS.append(copy(dust_sheets))
VARIABLES.append(copy(variables))
PLOTNAMES.append(copy(plotname))
ROWS.append(copy(rows))



'''main program'''
for i in range(len(PLOTNAMES)):
    print 'Producing plot %s' %PLOTNAMES[i]
    plot_dust(DUST_SHEETS[i], VARIABLES[i], PLOTNAMES[i], ROWS[i])


       
p.show()

