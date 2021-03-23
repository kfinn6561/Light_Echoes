'''
Created on Sep 23, 2015

@author: kieran
'''
import pylab as p
import numpy as np
from general_tools import pload,make_steps,pdump
import matplotlib.lines as mlines
import sys
from colour_tools import kelly_colors

import matplotlib as mpl
mpl.rcParams.update({'font.size': 24})
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']#Computer Modern Roman']
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.labelsize'] = 18.
mpl.rcParams['ytick.labelsize'] = 18.
mpl.rcParams['xtick.major.size']= 10.
mpl.rcParams['xtick.minor.size']= 5.
mpl.rcParams['ytick.major.size']= 10.
mpl.rcParams['ytick.minor.size']= 5.
mpl.rcParams['figure.max_open_warning']=False
mpl.rcParams['lines.linewidth']=2


extension='png'

data_folder='plot_data/'
plot_folder='plots/'
mag_base=100.**(1./5.)#base of the magnitude system

interesting_sne={'He':['Cas A','sn2011ei','sn1993J'],
                     'H_alpha':['Cas A','sn2008ax','sn1993J'],
                     'Ca_III':['Cas A','sn2003bg','sn2011ei'],
                     'Fe_II':['Cas A','sn2003bg','sn1996cb']}

maximum_sne={'He':['sn2011ei','sn1993J'],
                     'H_alpha':['sn2008ax','sn1993J'],
                     'Ca_III':[]}

line_names=['He','H_alpha','Ca_III']
les=['le3923','le2521','le2116']
line_plotnames={'He':'He I','H_alpha':r'H$\alpha$','Ca_III':'Ca II','Fe_II':'Fe 5169'}
le_plotnames={'le3923':'NE','le2521':'NW','le2116':'SW','unmodified':'w(t)=1'}

template_color='b'
#casa_color='#CD0000'
casa_color='k'
try:
    sn_colours=pload('sn_colors.pkl')
except IOError:
    sn_colours={'Cas A':casa_color,'Template':template_color}
kelly_colours_it=iter(kelly_colors)

sn_colours['sn2011ei_no33']='#FF6700'

selective=not ('all_sn' in sys.argv)
talk=('talk' in sys.argv)
flux=('flux' in sys.argv)
multifig=('multifig' in sys.argv)


all_plots=['features','windows','casa','bands_lines','individual_velocity','velocity','velocity_2',
           'velocity_3','velocity_4','velocity_5']
to_plot=all_plots
if 'individual' in sys.argv:
    to_plot.append('individual')
for item in sys.argv:#if any of the plot names appear in argv then only run those plots, else run all
    if item in all_plots:
        to_plot=sys.argv


    
    '''#####################################################################'''

if 'bands_lines' in to_plot:
    '''#####################################################################'''
    '''band lines plot'''
    print 'making the bands lines plot'
    p.figure(figsize=(10,9))
    band_colors={'B':'b','V':'g','r':'r','i':'m'}
    line_colors={'H_alpha':'b','He':'c','Fe_II':'g','Ca_III':'r'}
    data=pload(data_folder+'bands_lines.pkl')
    NN=0
    for datum in data['bands']:
        NN+=1
        band,x,y=datum
        p.plot(x,y,label=band,color=band_colors[band])
    
    for datum in data['lines']:
        NN+=1
        name,minl,maxl=datum
        p.axvline(minl,color=line_colors[name])
        p.axvline(maxl,color=line_colors[name])
        p.axvspan(minl,maxl,alpha=0.4,label=line_plotnames[name],color=line_colors[name])
    p.legend(loc='upper center', bbox_to_anchor=(0.5, 1),ncol=(NN+1)/2)
    p.xlabel('Wavelength (A)')
    p.ylabel('Flux')
    p.ylim(0,1.15)
    p.xlim(3050,9500)
    p.savefig(plot_folder+'bands_lines.pdf')
    '''#####################################################################'''


if 'windows' in to_plot:
    '''#####################################################################'''
    '''window function plot'''
    print 'making window function plot'
    colours=['k','r','b','y']
    data=pload(data_folder+'window_functions.pkl')
    sn93j_data=pload(data_folder+'sn1993J_window_functions.pkl')
    labels=['unmodified','le2116','le2521','le3923']
    temp=[le_plotnames[l] for l in labels]
    labels=temp
    alphas=[1,1,0.8,0.8]
    
    if talk:
        f_w,ax_w=p.subplots(2,1,sharex=True,figsize=(12,8))
        Nc=4
    else:
        f_w,ax_w=p.subplots(2,1,sharex=True,figsize=(7,7))
        Nc=2
    ts,windows=data['windows']
    
    p.tight_layout(h_pad=-1)
    for i in range(len(windows)):
        ax_w[0].plot(ts,windows[i],colours[i],lw=2)
    ax_w[0].set_ylim(0,1.1)
    ax_w[0].set_ylabel('w(t) Fractional Throughput')
    
    ts,light_curves=data['light_curves']
    sn93j_ts,sn93j_lc=sn93j_data['light_curves']
    for i in range(len(light_curves)):
        if flux:
            ax_w[1].plot(ts,mag_base**(-light_curves[i]),colours[i],lw=2)
        else:
            ax_w[1].plot(ts,light_curves[i],colours[i],lw=2,label=labels[i],alpha=alphas[i])
            ax_w[1].plot(sn93j_ts,sn93j_lc[i],colours[i],lw=2,ls='--')
        
    if flux:
        ax_w[1].set_ylabel('Flux')
        ax_w[1].set_ylim(0,1)
        ax_w[1].set_yticks([0.,0.2,0.4,0.6,0.8])
    else:
        ax_w[1].set_ylim(6,-1)
        ax_w[1].set_yticks([6,5,4,3,2,1,0])
        
        ax_w[1].set_ylabel('V Magnitude')
    
    ax_w[1].set_xlim(-15,95)
    ax_w[1].set_xlabel('Time (days)')
    
    ax_w[1].legend(loc='upper right', bbox_to_anchor=(1., 1.),ncol=Nc,fontsize=16)
    '''if not talk:
        b,t=ax_w[0].get_ylim()
        w=t-b
        ax_w[0].set_ylim(b,b+1.3*w)#increase the width so that the legend can be seen'''
    p.subplots_adjust(left=0.12,bottom=0.08)
    if flux:
        p.savefig(plot_folder+'windows_flux.png')
    else:
        p.savefig(plot_folder+'windows.pdf')
    '''#####################################################################'''



if 'individual' in to_plot:
    '''#####################################################################'''
    '''individual plot'''
    print 'making the individual plot'
    colours=['r','b','y','k']#last colour is for unmodified
    labels=['le2116','le2521','le3923','unmodified']
    temp=[le_plotnames[l] for l in labels]
    labels=temp
    data=pload(data_folder+'individual.pkl')
    for datum in data:
        sn_name,line_name,minl,maxl,lambdas,spectra=datum
        p.figure()
        for i in range(len(spectra)):
            x,y,z=spectra[i]
            p.plot(lambdas,x,colours[i],label=labels[i])
            p.plot(lambdas,y,colours[i]+'--')
            p.plot(lambdas,z,colours[i]+'--')
        p.xlim(minl-50,maxl+50)
        ax=p.gca()
        y_text=0.05
        ax.text(.5,y_text,'%s line for supernova %s' %(line_plotnames[line_name],sn_name),horizontalalignment='center',transform=ax.transAxes)
        p.xlabel('Wavelength (A)')
        p.ylabel('Flux')
        p.savefig(plot_folder+'spectra_%s_%s.pdf' %(sn_name,line_name))
    '''#####################################################################'''
 
 
mpl.rcParams['font.size']= 30#increase font size
mpl.rcParams['axes.labelsize']=30
mpl.rcParams['xtick.labelsize']=26
mpl.rcParams['ytick.labelsize']=26
       
if 'features' in to_plot:
    '''#####################################################################'''
    '''features plot'''
    print 'making the features plot'
    colours=['r','b','y','k']#last colour is for unmodified
    labels=['le2116','le2521','le3923','unmodified']
    temp=[le_plotnames[l] for l in labels]
    labels=temp
    data=pload(data_folder+'features.pkl')
    for datum in data:
        name,band,minl,maxl,lambdas,spectra=datum
        p.figure()
        for i in range(len(spectra)-int(talk)):
            x,y,z=spectra[i]
            p.plot(lambdas,x,colours[i],label=labels[i],lw=2)
            if i==len(spectra)-1:
                #p.plot(lambdas,y,colours[i]+'--')
                #p.plot(lambdas,z,colours[i]+'--')
                p.fill_between(lambdas,y,z,color=colours[i],alpha=0.4)
        p.xlim(minl-50,maxl+50)
        if name=='Ca_III':
            if talk:
                p.legend(loc='upper left',ncol=2,columnspacing=0.2,
                         handletextpad=0.1,borderpad=0.1,frameon=False)
            else:
                p.legend(loc='upper left')
            p.ylim(0.5,1.65)
        elif name=='H_alpha':
            a=mlines.Line2D([],[],color='k',label='Mean')
            b=mlines.Line2D([],[],color='k',ls='--',label=r'$\pm1\sigma$')
            #p.legend(handles=[a,b],loc='upper right')
            p.ylim(0.4,1.45)
        else:
            p.ylim(0.45,1.45)
            p.xlim(5320,6050)
        ax=p.gca()
        '''
        if name=='H_alpha':
            y_text=0.9
        else:
            y_text=0.05
        '''
        y_text=0.05
        if talk:
            ax.text(.05,y_text,line_plotnames[name],horizontalalignment='left',transform=ax.transAxes,fontsize=30)
        else:
            ax.text(.5,y_text,'%s line, calibrated with %s band' %(line_plotnames[name],band),horizontalalignment='center',transform=ax.transAxes)
        p.xlabel(r'Wavelength ($\AA$)')
        if name=='H_alpha' or (name=='He' and band=='r'):
            p.ylabel(' ')
            ax.yaxis.tick_right()
        else:
            p.ylabel(r'Scaled $f_\lambda$')
        p.tight_layout(h_pad=-1)
        p.savefig(plot_folder+'spectra_%s_%s.%s' %(name,band,extension))
    '''#####################################################################'''
                



if 'casa' in to_plot:
    '''#####################################################################'''
    '''casa plot'''
    print 'making the casa plot'
    
    le_color=template_color
    
    les,lines,data=pload(data_folder+'casa.pkl')
    
    
    if multifig:
        ax=[]
        figs=[]
        for i in range(len(lines)):
            ax.append([])
            figs.append([])
            for j in range(len(les)):
                ff,a=p.subplots(1,1)
                figs[i].append(ff)
                ax[i].append(a)
        ax=np.array(ax)
                
    else:
        if talk:
            f,ax=p.subplots(len(les),len(lines),sharex='col',sharey='row',figsize=(18.,12.))
        else:
            f,ax=p.subplots(len(les),len(lines),sharex='col',sharey='row',figsize=(17.,16.5))
        p.tight_layout(pad=3,h_pad=-1,w_pad=-1)
        ax=ax.T
    
    for axx in ax.flatten():
        axx.minorticks_on()
    ax2=np.ones_like(ax)
    for i in range(len(lines)):
        for j in range(len(les)):
            ax2[i][j]=ax[i][j].twiny()
            ax2[i][j].minorticks_on()
    
    
    for i in range(len(lines)):
        for j in range(len(les)):
            velo,le,casa=data[i][j]
            velo_labels,velo_positions=velo
            le_X,le_mid,le_min,le_max,le_x0,le_xmin,le_xmax=le
            casa_X,casa_mid,casa_min,casa_max,casa_x0,casa_xmin,casa_xmax=casa
            
            le_X=le_X/1000.
            le_x0/=1000.
            le_xmax/=1000.
            le_xmin/=1000.#rescale x axis
            casa_X=casa_X/1000.
            casa_x0/=1000.
            casa_xmax/=1000.
            casa_xmin/=1000.#rescale x axis
            
            velo_positions=velo_positions/1000.
            
            ax2[i][j].set_xticks(velo_positions)
            if j==0 or multifig:
                ax2[i][j].set_xticklabels([str(int(np.round(vl))) for vl in velo_labels])
                if i==1:
                    ax2[i][j].set_xlabel(r'Velocity ($10^3$ km/s)')
            else:
                ax2[i][j].set_xticklabels([])
            
            _,le_mid=make_steps(le_X,le_mid)
            _,le_min=make_steps(le_X,le_min)
            le_X,le_max=make_steps(le_X,le_max)
            
            _,casa_mid=make_steps(casa_X,casa_mid)
            _,casa_min=make_steps(casa_X,casa_min)
            casa_X,casa_max=make_steps(casa_X,casa_max)
            
            ax[i][j].plot(le_X,le_mid,color=le_color,label='Template')
            ax[i][j].fill_between(le_X,le_min,le_max,color=le_color,alpha=0.4)
            ax[i][j].axvline(le_x0,color=le_color)
            ax[i][j].axvspan(le_xmin,le_xmax,facecolor=le_color,alpha=0.4)
            
            ax[i][j].plot(casa_X,casa_mid,color=casa_color,label='Cas A')
            ax[i][j].fill_between(casa_X,casa_min,casa_max,color=casa_color,alpha=0.3)
            ax[i][j].axvline(casa_x0,color=casa_color)
            ax[i][j].axvspan(casa_xmin,casa_xmax,facecolor=casa_color,alpha=0.3)
            
            ax[i][j].set_xlim(le_X[0],le_X[-1]-0.01)
            ax[i][j].set_ylim(0,2)
            if j!=0:
                ax[i][j].set_yticks([0,0.5,1,1.5])
            #ax[i][j].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            ax2[i][j].set_xlim(ax[i][j].get_xlim())
        if not multifig:   
            ax[i][0].text(.05,0.87,line_plotnames[line_names[i]],horizontalalignment='left',transform=ax[i][0].transAxes, fontsize=28)
      
    for j in range(len(les)):
        if multifig:
            for i in range(len(ax)):
                ax[i][j].set_ylabel(r'Scaled $f_\lambda$')
        else:                
            ax[0][j].set_ylabel(r'Scaled $f_\lambda$')
            ax[0][j].text(.05,0.1,le_plotnames[les[j]],horizontalalignment='left',transform=ax[0][j].transAxes, fontsize=28)
    if multifig:
        for i in range(len(ax)):
            for j in range(len(les)):
                ax[i][j].legend(loc='upper left',bbox_to_anchor=(0.01,0.96))                  
    else:
        ax[-1][-1].legend(loc='upper left',bbox_to_anchor=(0.01,0.96),fontsize=24)
    ax[1][-1].set_xlabel(r'Wavelength ($10^3\AA$)')
    if multifig:
        for i in range(len(ax)):
            for j in range(len(les)):
                ax[i][j].set_xlabel(r'Wavelength ($10^3\AA$)')
    if multifig:
        for i in range(len(figs)):
            for j in range(len(les)):
                ax[i][j].text(.05,0.1,line_plotnames[line_names[i]],horizontalalignment='left',transform=ax[i][j].transAxes, fontsize=24)
                ax[i][j].text(.95,0.1,les[j],horizontalalignment='right',transform=ax[i][j].transAxes, fontsize=24)
                ax[i][j].set_ylim(0.5,1.5)
                figs[i][j].savefig(plot_folder+'casa_%s_%s.png' %(line_names[i],les[j]))
    p.savefig(plot_folder+'casa.%s' %extension)
    '''#####################################################################'''




if 'individual_velocity' in to_plot:
    '''#####################################################################'''
    '''individual velocity plot'''
    print 'making the individual velocity plot'
    
    
    les,lines,data=pload(data_folder+'individual_velocity.pkl')
    _,_,casa_data=pload(data_folder+'casa.pkl')
    xlims={'He':(5.56,5.84),'H_alpha':(6.12,6.44),'Ca_III':(7.8,8.8)}
    
    if talk:
        f,ax=p.subplots(len(les),len(lines),sharex='col',sharey='row',figsize=(18,12))
    else:
        f,ax=p.subplots(len(les),len(lines),sharex='col',sharey='row',figsize=(17,21))
    p.tight_layout(pad=3,h_pad=-1,w_pad=-1)
    for axx in ax.flatten():
        axx.minorticks_on()
    ax=ax.T
    
    ax2=np.ones_like(ax)
    for i in range(len(lines)):
        for j in range(len(les)):
            ax2[i][j]=ax[i][j].twiny()
    
    
    handles={}
    for i in range(len(lines)):
        for j in range(len(les)):
            velo,sne=data[i][j]
            _,template,casa=casa_data[i][j]
            sne=[['Cas A']+casa,['Template']+template]+sne
            
            '''
            temp=[]
            for sn in sne:
                if sn[0] in ['sn2011ei','sn2011ei_no33']:
                    temp.append(sn)
            sne=temp
            '''
                
            velo_labels,velo_positions=velo
            
            ax2[i][j].set_xticks(velo_positions)
            if j==0:
                ax2[i][j].set_xticklabels([str(int(np.round(vl))) for vl in velo_labels])
                ax2[i][j].set_xlabel(r'Velocity ($10^3$ km/s)')
            else:
                ax2[i][j].set_xticklabels([])
            
            for sn in sne:
                name,le_X,le_mid,le_min,le_max,le_x0,le_xmin,le_xmax=sn
                le_X=le_X/1000.
                le_x0/=1000.
                le_xmax/=1000.
                le_xmin/=1000.#rescale x axis
                if selective and (name not in interesting_sne[lines[i]]):
                    continue
                try:
                    le_color=sn_colours[name]
                except KeyError:
                    sn_colours[name]=next(kelly_colours_it)
                    le_color=sn_colours[name]
                handles[name]=ax[i][j].plot(le_X,le_mid,color=le_color,label=name)[0]
                if name !='Cas A':
                    ax[i][j].fill_between(le_X,le_min,le_max,color=le_color,alpha=0.4)
                elif selective:
                    ax[i][j].fill_between(le_X,le_min,le_max,color=le_color,alpha=0.2)
                ax[i][j].axvline(le_x0,color=le_color)
                if name=='Cas A':
                    ax[i][j].axvspan(le_xmin,le_xmax,facecolor=le_color,alpha=0.2)
                else:
                    ax[i][j].axvspan(le_xmin,le_xmax,facecolor=le_color,alpha=0.4)
            
            lo,hi=xlims[line_names[i]]
            ax[i][j].set_xlim(lo,hi)
            ax[i][j].set_ylim(0,2)
            if j!=0:
                ax[i][j].set_yticks([0,0.5,1,1.5])
            #ax[i][j].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            ax2[i][j].set_xlim(ax[i][j].get_xlim())
           
        ax[i][0].text(.05,0.1,line_plotnames[line_names[i]],horizontalalignment='left',transform=ax[i][0].transAxes, fontsize=24)
    hnames=handles.keys()
    try:
        hnames.remove('Cas A')
        hnames=['Cas A']+hnames
    except ValueError:
        pass
    hvals=[handles[hn] for hn in hnames]  
    for j in range(len(les)):
        ax[0][j].set_ylabel(r'Scaled $f_\lambda$')
        ax[0][j].text(.05,0.9,le_plotnames[les[j]],horizontalalignment='left',transform=ax[0][j].transAxes, fontsize=24)
    ax[-1][0].legend(handles=hvals,loc='upper left',bbox_to_anchor=(0.01,0.96))
    for i in range(len(lines)):
        ax[i][-1].set_xlabel(r'Wavelength ($10^3\AA$)')
    if selective:
        p.savefig(plot_folder+'individual_velocity_selective.pdf')
    else:
        p.savefig(plot_folder+'individual_velocity.pdf')
    '''#####################################################################'''


   
if 'velocity' in to_plot:
    '''#####################################################################'''
    '''velocity plot'''
    print 'making the velocity plot'
    data=pload(data_folder+'velocity_list.pkl')
    
    lines=line_names[:-1]
    markers=['o','*','v','+']
    le_markers={les[i]:markers[i] for i in range(len(les))}
    #sne=data[lines[0]][les[0]].keys()#if the order is funky then may want to change this
    sne=['sn2011ei','sn2003bg','sn2008ax','Cas A','sn1996cb','sn1993J','sn2011dh']
    dline=35
    dle=2.5
    dsn=10
    p.figure(figsize=(12,7))
    ticks=[]
    labels=[]
    ms=10
    
    for i in range(len(lines)):
        line=lines[i]
        temp_sne=[]
        for sn in sne:#this is probably a convoluted why to get the sne in the right order, but it should do the trick
            if sn in interesting_sne[line]:
                temp_sne.append(sn)
        for k in range(len(temp_sne)):
            sn=temp_sne[k]
            try:
                sn_color=sn_colours[sn]
            except KeyError:
                sn_colours[sn]=next(kelly_colours_it)
                sn_color=sn_colours[sn]
            ticks.append(i*dline+k*dsn+dle)
            labels.append(sn)
            for j in range(len(les)):
                le=les[j]
                x=i*dline+j*dle+k*dsn
                y,ymin,ymax=data[line][le][sn]
                p.errorbar(x,y,yerr=[[y-ymax],[ymin-y]],ls='None',marker=le_markers[le],markersize=ms,color=sn_color)
    
    for le in les:
        p.plot([],[],ls='None',marker=le_markers[le],markersize=ms,color='k',label=le_plotnames[le])#dummy lines to provide labels for the legend
    p.legend(loc='lower center', bbox_to_anchor=(0.5, 0),ncol=3,fontsize=18,numpoints=1)
    p.xlim(-dle,(len(lines)-1)*dline+len(interesting_sne[lines[-1]])*dsn)
    p.xticks(ticks,rotation=45)
    ax=p.gca()
    ax.set_xticklabels(labels)
    for i in range(len(lines)):
        ax.text(i*dline,-14.5,line_plotnames[lines[i]],horizontalalignment='left')
    p.ylim(-4.1,-15.5)
    p.ylabel(r'Velocity ($10^3$ km/s)')
    p.subplots_adjust(left=0.08,bottom=0.2,right=0.95)
    p.savefig(plot_folder+'velocity.%s' %extension)
                
    
    '''#####################################################################'''
    
    
if 'velocity_3' in to_plot:
    '''#####################################################################'''
    '''velocity 3 plot'''
    print 'making the velocity 3 plot'
    data=pload(data_folder+'velocity_list.pkl')
    
    
    lines=line_names#[:-1]
    #sne=data[lines[0]][les[0]].keys()#if the order is funky then may want to change this
    sne=['Cas A','Template','sn1993J','sn2003bg','sn1996cb','sn2008ax','sn2011dh','sn2011ei']
    filled=['Cas A','Template','sn2003bg','sn1993J']
    markers={sn:'o' for sn in sne}
    markers['Cas A']='s'
    markers['Template']='d'

    dle=28
    dsn=2.5
    f,ax=p.subplots(len(lines),2,sharex=False,figsize=(16,8),gridspec_kw={'width_ratios':[4,1]})
    ax_v=ax.T[0]
    ax_d=ax.T[1]
    f.delaxes(ax_d[-1])
    ax_d=ax_d[:-1]

    
    p.tight_layout(h_pad=-1)
    ticks=[]
    labels=[]
    ms=10
    
    for i in range(len(lines)):
        line=lines[i]
        for j in range(len(les)):
            le=les[j]
            temp_sne=sorted(data[line][le].keys(),key=lambda x:data[line][le][x][0])
            if i==0:
                ticks.append(len(temp_sne)/2.*dsn+j*dle)
                labels.append(le_plotnames[le])
            d_mins=[]
            d_maxs=[]
            for k in range(len(temp_sne)):
                sn=temp_sne[k]
                try:
                    sn_color=sn_colours[sn]
                except KeyError:
                    sn_colours[sn]=next(kelly_colours_it)
                    sn_color=sn_colours[sn]
                x=j*dle+k*dsn
                y,ymin,ymax=data[line][le][sn]
                
                if sn in filled:
                    mfc=sn_color
                else:
                    mfc='none'
                ax_v[i].errorbar(x,y,yerr=[[y-ymax],[ymin-y]],ls='None',marker=markers[sn],markersize=ms,
                                 mec=sn_color,mew=1.5,mfc=mfc,color=sn_color)
                if sn in ['sn1993J','sn2003bg']:
                    ax_v[i].errorbar(x,y,yerr=[[y-ymax],[ymin-y]],ls='None',mew=1.5,color='k')
                if sn in maximum_sne[line]:
                    casa_mid,casa_min,casa_max=np.array(data[line][le]['Cas A'])
                    casa_err_n=casa_min-casa_mid
                    casa_err_p=casa_mid-casa_max
                    sn_err_n=ymin-y
                    sn_err_p=y-ymax
                    
                    err_n=np.sqrt(casa_err_p**2+sn_err_n**2)
                    err_p=np.sqrt(casa_err_n**2+sn_err_p**2)
                    d_mid=casa_mid-y
                    
                    d_mins.append(d_mid-2*err_n)
                    d_maxs.append(d_mid+2*err_p)
                    
                    ax_d[i].errorbar(dle*j,d_mid,yerr=[[err_n],[err_p]],ls='None',marker=markers[sn],
                                     mec=sn_color,mew=1.5,mfc=mfc,markersize=ms,color=sn_color)
            if line!='Ca_III':
                x1=j*dle-dle/4.
                x2=j*dle+dle/4.
                y1=max(d_maxs)
                y2=min(d_mins)
                ax_d[i].plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],color='k',lw=1)
                    
    
    for sn in sne:
        if sn in filled:
            mec='none'
            mew=0
            mfc=sn_colours[sn]
        else:
            mec=sn_colours[sn]
            mew=1.5
            mfc='none'
        p.plot([],[],ls='None',marker=markers[sn],markersize=ms,
               mec=mec,mew=mew,mfc=mfc,color=sn_colours[sn],label=sn)#dummy lines to provide labels for the legend
    ax_v[-1].legend(loc='lower center', bbox_to_anchor=(0.57, 0),ncol=4,fontsize=23,numpoints=1,columnspacing=0.8,handletextpad=0.1)
    
    for i in range(len(ax_v)):
        ax_v[i].set_xticks(ticks)
        ax_v[i].set_xticklabels(['']*len(ticks))
        ax_v[i].set_xlim(-dsn,(len(les)-1)*dle+len(sne)*dsn)
    ax_v[-1].set_xticklabels(labels,fontsize=28)
    
    ax_v[1].set_ylabel(r'Velocity ($10^3$ km/s)')
    ax_v[0].set_ylim(-6.1,-13.1)
    ax_v[1].set_ylim(-8.1,-14.5)
    #ax_v[2].set_ylim(-4.1,-12.5)
    ax_v[2].set_ylim(-2.4,-12.5)
    
    for i in range(len(lines)):
        ax_v[i].text(0.05,0.87,line_plotnames[lines[i]],horizontalalignment='left',transform=ax_v[i].transAxes,fontsize=28)    
        
    for axis in ax_d:
        axis.yaxis.tick_right()
        axis.yaxis.set_ticks_position('both')
        axis.axhline(y=0.,color='k')
        axis.set_ylim(-8,8)
        axis.set_yticks([-4,0,4])
        axis.set_xticks(dle*np.arange(len(les)))
        axis.set_xticklabels(['']*len(les))
        axis.set_xlim(-dle/2.,dle*(len(les)-0.5))
    
    ax_d[-1].set_xticklabels(labels,fontsize=28)    
    ax_d[1].set_ylabel(r'Velocity difference ($10^3$ km/s)',rotation=270,labelpad=20)
    ax_d[1].yaxis.set_label_position('right')
    ax_d[1].yaxis.set_label_coords(1.3,0.5)
    p.subplots_adjust(left=0.08,bottom=0.08,right=0.9,wspace=0.,hspace=0.)
    
    p.savefig(plot_folder+'velocity_3.pdf')
    
'''###############################################################'''
    
if 'velocity_2' in to_plot:
    '''#####################################################################'''
    '''velocity 2 plot'''
    print 'making the velocity 2 plot'
    data=pload(data_folder+'velocity_list.pkl')
    
    
    lines=line_names#[:-1]
    marker='o'
    #sne=data[lines[0]][les[0]].keys()#if the order is funky then may want to change this
    sne=['Cas A','Template','sn2011ei','sn2003bg','sn2008ax','sn1996cb','sn1993J','sn2011dh']
    dle=28
    dsn=2.5
    f_v,ax_v=p.subplots(len(lines),1,sharex=True,figsize=(12,8))
    p.tight_layout(h_pad=-1)
    ticks=[]
    labels=[]
    ms=10
    
    for i in range(len(lines)):
        line=lines[i]
        for j in range(len(les)):
            le=les[j]
            temp_sne=sorted(data[line][le].keys(),key=lambda x:data[line][le][x][0])
            if i==0:
                ticks.append(len(temp_sne)/2.*dsn+j*dle)
                labels.append(le_plotnames[le])
            for k in range(len(temp_sne)):
                sn=temp_sne[k]
                try:
                    sn_color=sn_colours[sn]
                except KeyError:
                    sn_colours[sn]=next(kelly_colours_it)
                    sn_color=sn_colours[sn]
                x=j*dle+k*dsn
                y,ymin,ymax=data[line][le][sn]
                ax_v[i].errorbar(x,y,yerr=[[y-ymax],[ymin-y]],ls='None',marker=marker,markersize=ms,color=sn_color)
    
    for sn in sne:
        p.plot([],[],ls='None',marker=marker,markersize=ms,color=sn_colours[sn],label=sn)#dummy lines to provide labels for the legend
    p.legend(loc='lower center', bbox_to_anchor=(0.57, 0),ncol=4,fontsize=18,numpoints=1)
    p.xlim(-dsn,(len(les)-1)*dle+len(sne)*dsn)
    ax_v[-1].set_xticks(ticks)
    ax_v[-1].set_xticklabels(labels)
    
    ax_v[1].set_ylabel(r'Velocity ($10^3$ km/s)')
    ax_v[0].set_ylim(-6.1,-13.1)
    ax_v[1].set_ylim(-8.1,-14.5)
    ax_v[2].set_ylim(-4.1,-12.5)
    if talk:
        p.subplots_adjust(left=0.11,bottom=0.08,right=0.95)
    else:
        p.subplots_adjust(left=0.08,bottom=0.08,right=0.95)
    for i in range(len(lines)):
        if talk:
            ax_v[i].text(0.05,0.85,line_plotnames[lines[i]],horizontalalignment='left',transform=ax_v[i].transAxes)
        else:
            ax_v[i].text(0.05,0.9,line_plotnames[lines[i]],horizontalalignment='left',transform=ax_v[i].transAxes)
    p.savefig(plot_folder+'velocity_2.%s' %extension)
    
'''###############################################################'''
    
if 'velocity_4' in to_plot:
    '''#####################################################################'''
    '''velocity 4 plot'''
    print 'making the velocity 4 plot'
    
    data=pload(data_folder+'velocity_list.pkl')

    lines=line_names[:-1]
    markers=['o','*','v','+']
    le_markers={les[i]:markers[i] for i in range(len(les))}
    #sne=data[lines[0]][les[0]].keys()#if the order is funky then may want to change this
    sne=['sn2011ei','sn2003bg','sn2008ax','Cas A','sn1996cb','sn1993J','sn2011dh']
    dline=5.
    dle=1.
    dsn=0
    #p.figure(figsize=(7.5,8))
    p.figure(figsize=(7.5,10))
    ticks=[]
    labels=[]
    ms=15
    
    for i in range(len(lines)):
        line=lines[i]
        temp_sne=[]
        for sn in sne:#this is probably a convoluted why to get the sne in the right order, but it should do the trick
            if sn in interesting_sne[line]:
                temp_sne.append(sn)
        for k in range(len(temp_sne)):
            sn=temp_sne[k]
            try:
                sn_color=sn_colours[sn]
            except KeyError:
                sn_colours[sn]=next(kelly_colours_it)
                sn_color=sn_colours[sn]
            ticks.append(i*dline+k*dsn+dle)
            labels.append(sn)
            for j in range(len(les)):
                le=les[j]
                x=i*dline+j*dle+k*dsn
                y,ymin,ymax=data[line][le][sn]
                p.errorbar(x,y,yerr=[[y-ymax],[ymin-y]],ls='None',marker=le_markers[le],capsize=5,capthick=2,markersize=ms,color=sn_color,alpha=0.5)
                p.plot(x,y,ls='None',marker=le_markers[le],markersize=ms,color=sn_color)
                
    
    ax=p.gca()
    
    leg_le_x1=1.
    leg_le_dx=2.2
    leg_le_y=-4.7
    leg_dxname=0.4
    leg_sn_x1=1.1
    leg_sn_x1=leg_le_x1
    leg_sn_dx=3.5
    leg_dy=0.7
    leg_bounding=0.3
    leg_rbounding=0.35
    
    for i in range(len(les)):
        le=les[i]
        p.plot([leg_le_x1+leg_le_dx*i],[leg_le_y],ls='None',marker=le_markers[le],markersize=10,color='k')
        ax.text(leg_le_x1+leg_dxname+leg_le_dx*i,leg_le_y,le_plotnames[le],fontsize=28,ha='left',va='center')
        
    
    plot_sne=['Cas A','sn1993J','sn2008ax','sn2011ei']
    for i in range(len(plot_sne)):
        sn=plot_sne[i]
        yi=2-i%2
        xi=i/2
        p.plot([leg_sn_x1+xi*leg_sn_dx],[leg_le_y-yi*leg_dy],ls='None',marker='o',markersize=10,color=sn_colours[sn])
        ax.text(leg_sn_x1+xi*leg_sn_dx+leg_dxname,leg_le_y-yi*leg_dy,sn,fontsize=28,ha='left',va='center')
    
    xbox=[leg_le_x1-leg_bounding,
          leg_sn_x1+2*leg_sn_dx-leg_rbounding,
          leg_sn_x1+2*leg_sn_dx-leg_rbounding,
          leg_sn_x1-leg_bounding,
          leg_sn_x1-leg_bounding,
          leg_le_x1-leg_bounding,
          leg_le_x1-leg_bounding]
    ybox=[leg_le_y+0.5*leg_dy,
          leg_le_y+0.5*leg_dy,
          leg_le_y-2.5*leg_dy,
          leg_le_y-2.5*leg_dy,
          leg_le_y-0.5*leg_dy,
          leg_le_y-0.5*leg_dy,
          leg_le_y+0.5*leg_dy]
    
    p.plot(xbox,ybox,lw=1,color='k')
    
    
    p.xlim(-dle,dline+3*dle)
    p.xticks([])
    
    #ax.set_xticklabels(labels)
    for i in range(len(lines)):
        ax.text(i*dline,-14.5,line_plotnames[lines[i]],horizontalalignment='left',fontsize=28)
    #p.ylim(-5.8,-15.5)
    p.ylim(-4.2,-15.5)
    p.ylim(-3.9,-15.5)
    
    p.ylabel(r'Velocity ($10^3$ km/s)')
    p.subplots_adjust(left=0.13,top=0.95,bottom=0.05,right=0.95)
    p.savefig(plot_folder +'velocity_4.pdf')
    
    
if 'velocity_5' in to_plot:
    '''#####################################################################'''
    '''velocity 5 plot'''
    print 'making the velocity 5 plot'
    
    data=pload(data_folder+'velocity_list.pkl')

    lines=line_names[:-1]
    markers=['o','*','v','+']
    le_markers={les[i]:markers[i] for i in range(len(les))}
    #sne=data[lines[0]][les[0]].keys()#if the order is funky then may want to change this
    sne=['sn2011ei','sn2003bg','sn2008ax','Cas A','sn1996cb','sn1993J','sn2011dh']
    dline=5.
    dle=1.
    p.figure(figsize=(7.5,8))
    ms=10
    
    sne_handles={}
    for i in range(len(lines)):
        line=lines[i]
        for j in range(len(les)):
            le=les[j]
            x=i*dline+j*dle
            xbox=[x-0.3*dle,x+0.3*dle,x+0.3*dle,x-0.3*dle,x-0.3*dle]
            ybox=[]
            for sn in interesting_sne[line]:
                try:
                    sn_color=sn_colours[sn]
                except KeyError:
                    sn_colours[sn]=next(kelly_colours_it)
                    sn_color=sn_colours[sn]
                y,ymin,ymax=data[line][le][sn]
                    
                if sn=='Cas A':
                    p.errorbar(x,y,yerr=[[y-ymax],[ymin-y]],ls='None',marker=le_markers[le],markersize=ms,color=sn_color,alpha=0.5)
                    p.plot(x,y,ls='None',marker=le_markers[le],markersize=ms,color=sn_color)
                else:
                    p.errorbar(x,y,yerr=[[y-ymax],[ymin-y]],ls='None',marker='None',color='k')
                    ybox+=[y,y]
            ybox.append(ybox[0])
            p.plot(xbox,ybox,'k')
    
    le_handles=[]
    for le in les:
        le_handles+=p.plot([],[],ls='None',marker=le_markers[le],markersize=ms,color='k',label=le_plotnames[le])#dummy lines to provide labels for the legend
    p.legend(handles=le_handles,loc='lower right', bbox_to_anchor=(0.99, 0),ncol=3,fontsize=28,numpoints=1,columnspacing=0.8,handletextpad=0.1)
    
    p.xlim(-dle,dline+3*dle)
    p.xticks([])
    ax=p.gca()
    #ax.set_xticklabels(labels)
    for i in range(len(lines)):
        ax.text(i*dline,-14.5,line_plotnames[lines[i]],horizontalalignment='left',fontsize=28)
    p.ylim(-5.8,-15.5)
        
    p.ylabel(r'Velocity ($10^3$ km/s)')
    p.subplots_adjust(left=0.13,top=0.95,bottom=0.05,right=0.95)
    p.savefig(plot_folder +'velocity_5.pdf')


pdump(sn_colours,'sn_colors.pkl')

p.show()



