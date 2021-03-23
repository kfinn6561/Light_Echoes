import numpy as np
from scipy.io.idl import readsav
import matplotlib.pylab as pl
import pylabsetup
import sys
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import readcol 
import time
import glob
from matplotlib.backends.backend_pdf import PdfPages
import pickle as pkl
from colour_tools import kelly_colors

sne=['sn1998fa','sn2000H','sn2004ff','sn2006T','sn2006el','sn2008bo','sn2009mg','sn2011fu',
     'sn1996cb','sn2008ax','sn2011dh','sn2011ei','sn1993J','sn2003bg']
alphabetical_sne=['sn1998fa','sn2000H','sn2004ff','sn2006T','sn2006el','sn2008bo','sn2009mg','sn2011fu','sn1993J',
     'sn1996cb','sn2003bg','sn2008ax','sn2011dh','sn2011ei']

sne_handles={}

data_fol='yuqian_plot_data/'
colours=iter(kelly_colors[6:])
sne_colours={}


def plot_indi(vel_list,vel_dir,xr,yr,save_plot, annot = ''):
    # Kieran's color    
    sn_c = pkl.load(open("sn_colors.pkl",'rb'))
    key = sn_c.keys()

    fig,ax = pl.subplots(figsize=(7.5,7.5))
    pl.xlim(xr[0], xr[1])
    pl.ylim(yr[0], yr[1]) # in unit of 1000 km/s
        
    filename=open(data_fol+vel_list,'r').read()
    sn_name_list=filename.split('\n')
    symbol=['-<','->','-^','-v','-*','-d','-s','-p', '-h']
    j = 0
    sne_files={sn.split('_')[0]:sn for sn in sn_name_list}
    
    for i, sn in enumerate(sne):
        sn_name=sne_files[sn] 
        if any(sn_name.split('_')[0] in s for s in key): # to be consistent with color in Kieran's paper
            res = [x for x in key if sn_name.split('_')[0] in x]
            spec, phase, vel,velerr=readcol.readcol(vel_dir+sn_name,twod=False)
            MFC='none' if sn not in ['sn1993J', 'sn2003bg'] else sn_c[res[0]]
            sne_handles[sn]=ax.errorbar(phase, vel/1000, yerr=[velerr/1000, velerr/1000],capthick=2,fmt='-o',ms=6.5,
            label=sn_name.split('_')[0],color=sn_c[res[0]],mec=sn_c[res[0]],mfc=MFC,mew=1.5)
        else:
            spec, phase, vel,velerr=readcol.readcol(vel_dir+sn_name,twod=False)
            try:
                c=sne_colours[sn]
            except KeyError:
                c=next(colours)
                sne_colours[sn]=c
            sne_handles[sn]=ax.errorbar(phase, vel/1000, yerr=[velerr/1000, velerr/1000],capthick=2,fmt=symbol[j%9],mew=0,ms=8,
            label=sn_name.split('_')[0],color='gray')
            j = j+1
    pl.text((xr[1]-xr[0])*0.1+xr[0],(yr[1]-yr[0])*0.9+yr[0],annot,fontsize=20)
    pl.xlabel("Phase since V-band maximum (days)",fontsize=20)
    pl.ylabel("Absorption velocity ($10^3$ km s$^{-1}$)",fontsize=20)    
    minorLocatory   = MultipleLocator(1000)
    minorLocatorx   = MultipleLocator(10)
    ax.xaxis.set_minor_locator(minorLocatorx)
    ax.yaxis.set_minor_locator(minorLocatory)
    pl.legend(handles=[sne_handles[sn] for sn in alphabetical_sne],fontsize=15,mode="expand",loc=3,ncol=2,bbox_to_anchor=(0.3, .6, 0.7, 1))    
    pl.subplots_adjust(left=0.15)
    pl.savefig(save_plot)
    #pl.close()
    
plot_indi("inputIIb_HeI5875",data_fol,[-25,130],[-3,-17],'plots/IIb_HeI5876_vabs.pdf',annot = 'He I')
plot_indi("inputIIb_Halpha",data_fol,[-25,70],[-8,-25],'plots/IIb_Halpha_vabs.pdf', annot =r'H$\alpha$')
pl.show()