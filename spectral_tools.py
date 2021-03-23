'''
Created on May 27, 2015

@author: kieran
'''
#from integratespec_cfa3 import intspecclass
import numpy as np
import pylab as p
import os
from scipy.interpolate import interp1d
from scipy.io.idl import readsav
from colour_tools import linear_gradient
from general_tools import progress_bar,pload,pdump,reduce_range,find_nearest,mat_x_vec
from functions import  integrate,integrate_data,log_n,fit_spline,plot_spline,newton,gaussian
from copy import copy
from numpy.linalg.linalg import LinAlgError
    
mag_base=100.**(1./5.)#base of the magnitude system
SAV_FOLDER='meanspec/'
VERBOSE=False

zeropoints={}
f=open('band_pass_filters/zeropoints.txt','r')
for line in f:
    a,b=line.split()
    zeropoints[a]=log_n(float(b),mag_base)
f.close()
zeropoints['B-V']=zeropoints['B']-zeropoints['V']
zeropoints['V-r']=zeropoints['V']-zeropoints['r']
zeropoints['r-i']=zeropoints['r']-zeropoints['i']

def extrapolate(x,y,x_new,Npoints=500, xmin=50):#extends a list of x,y values to x_new by linear extrapolation after xmin
    x=np.array(x)
    y=np.array(y)
    idx=x>xmin
    if np.sum(idx)<=2:
        x1,x2=x[-2:]
        y1,y2=y[-2:]
        m=(y1-y2)/(x1-x2)
        c=y1-m*x1
    else:
        m,c=np.polyfit(x[idx],y[idx],1)
        
    
    x_new=np.linspace(x[-1],x_new,Npoints)[1:]
    y_new=m*x_new+c
    return(np.append(x,x_new),np.append(y,y_new))


class Band():
    def __init__(self,f,minl,maxl):
        self.f=f
        self.minl=minl
        self.maxl=maxl
    def __call__(self,x):
        return self.f(x)


class window_function():
    #a class which is essentially just a function but also has information about the support of the function
    def __init__(self,f,support):
        self.f=f
        self.support=support
    def __call__(self,x):
        return self.f(x)

'''
class spec_plotter(intspecclass):
    def extract_data(self):
        #extracts the data and returns it in a more usable format for plotting
        try:
            return [self.phase,self.wlength,self.flux] #check to see if it has already been run
        except AttributeError:
            phase_keys=np.sort(self.data.keys())
            phase=[]
            wlength=[]
            flux=[]
            for i in phase_keys:
                phase.append(self.data[i]['restphase'])
                data_keys=np.sort(self.data[i]['spec'].data.keys())
                temp_wlength=[]
                temp_flux=[]       
                for j in data_keys:
                    temp_wlength.append(self.data[i]['spec'].data[j]['lambda'])
                    temp_flux.append(self.data[i]['spec'].data[j]['flux'])
                wlength.append(np.array(temp_wlength))
                flux.append(np.array(temp_flux))
            self.phase=np.array(phase)
            self.wlength=np.array(wlength)
            self.flux=np.array(flux)
            return [self.phase,self.wlength,self.flux]
    
    
    def timeseries(self,offset=0.05*1e-12,lines=[10500.,11000.]):
        phase,wlength,flux=self.extract_data()
        plot_timeseries(phase, wlength, flux, offset, lines)
        
    def get_lightcurve(self,band='V'):
        try:
            return [self.lc_phases,self.lc_flux]#check if already run
        except AttributeError:
            if self.lc.data=={}:#light curve data has not been loaded
                self.loadlc_MLCS()
            phase_keys=np.sort(self.lc.data.keys())
            phases=[]
            fluxes=[]
            for i in phase_keys:
                phases.append(self.lc.data[i]['restphase'])
                fluxes.append(self.lc.data[i][band])
            self.lc_phases=np.array(phases)
            self.lc_flux=np.array(fluxes)
            lc_flux=mag_base**(-self.lc_flux)#lc_flux is in magnitudes
            self.lightcurve=interp1d(self.lc_phases,lc_flux)            
            return [self.lc_phases,self.lc_flux]
        
    def get_photospectra(self,band='V'):#combine the spectrum and photometry
        try:
            return self.photoflux
        except AttributeError:
            phase,wlength,flux=self.extract_data()
            lc_phase,lc_flux=self.get_lightcurve(band)
            self.photoflux,_=calc_photospectra(phase,wlength,flux,lc_phase,lc_flux,band)
            return self.photoflux
'''      
 
class Template():
    def __init__(self):
        self.photometry_data='colormeans_corrected.pkl'
        self.light_curve_data='light_curves.dat'
        self.extrapolate_lightcurve=False
        
    def get_lightcurve(self,band='V'):
        try:
            return [self.lc_phases[band],self.lc_flux[band]]#check if already run
        except AttributeError:
            lc_phases,lc_flux=pload(self.light_curve_data)
            self.t1=lc_phases
            self.t2=lc_flux
            self.lightcurve={}
            self.lightcurve_mags={}
            self.lc_phases={}
            self.lc_flux={}
            for b in lc_flux.keys():
                indices=~np.isnan(lc_flux[b])
                self.lc_phases[b]=lc_phases[indices]
                self.lightcurve_mags[b]=interp1d(self.lc_phases[b],copy(lc_flux[b][indices]),bounds_error=False,fill_value=np.inf)
                if self.extrapolate_lightcurve:
                    self.lc_phases[b],tflux=extrapolate(self.lc_phases[b],lc_flux[b][indices],self.extrapolate_lightcurve)
                else:
                    tflux=copy(lc_flux[b][indices])
                self.lc_flux[b]=mag_base**(-tflux)#lc flux is in magnitudes
                self.lightcurve[b]=interp1d(self.lc_phases[b],tflux,bounds_error=False,fill_value=0)#This will discard any spectra where ther is no light-curve data
                self.lightcurve_mags[b]=interp1d(self.lc_phases[b],tflux,bounds_error=False,fill_value=np.inf)
            return [self.lc_phases[band],self.lc_flux[band]]
        
    def get_photometry(self):
        try:
            return self.photometry
        except AttributeError:
            out={}
            photo=pload(self.photometry_data)
            for band in ['B-V','V-r','r-i']:
                out[band]=get_photo(photo[band],band)
            self.photometry=out
            return self.photometry
        
    def get_photo_data(self,fname):
        #pkl_name='photometry/%s_photo.pkl' %self.sn_name
        pkl_name=fname.split('.')[0]+'.pkl'
        self.light_curve_data=pkl_name
        if os.path.exists(pkl_name):
            return #pickle file already created            
        f=open(fname,'r')
        lines=f.readlines()
        bands=lines[0].split(',')[1:]
        bands=[b.strip() for b in bands]
        phase=[]
        flux={b:[] for b in bands}
        for line in lines[1:]:
            data=line.split(',')
            phase.append(float(data[0]))
            for i in range(1,len(data)):
                try:
                    flux[bands[i-1]].append(float(data[i].split()[0]))#must split because second entry is error
                except (ValueError, IndexError):
                    flux[bands[i-1]].append(np.nan)#store nan for any date that doesn't have observation.
        phase=np.array(phase)
        for b in bands:
            flux[b]=np.array(flux[b])
        pdump([phase,flux],pkl_name)
        return
        
    def get_photospectra(self,band='V'):#combine the spectrum and photometry
        try:
            return self.photoflux[band]
        except AttributeError:
            self.photoflux={}
            self.photomin={}
            self.photomax={}
            self.photonorms={}
        except KeyError:
            pass
        phase,wlength,flux=self.extract_data()
        lc_phase,lc_flux=self.get_lightcurve(band=band)
        photomin=self.error_n
        photomax=self.error_p
        if self.continuum=='cubic':
            out=[]
            print '\nadding continuum'
            for i in range(len(phase)):
                progress_bar(i, len(phase))
                continuum=self.calc_continuum(wlength[i],flux[i],phase[i])
                out.append((flux[i]+1.)*continuum)
                photomin[i]*=continuum
                photomax[i]*=continuum
            print '\n'
        elif self.continuum=='one':
            out=flux+1
        else:
            out=flux
        self.photoflux[band],norms=calc_photospectra(phase,wlength,out,lc_phase,lc_flux,band)
        self.photonorms[band]=norms
        self.photomin[band]=np.array([self.photoflux[band][i]-photomin[i]*norms[i] for i in range(len(norms))])
        self.photomax[band]=np.array([self.photoflux[band][i]+photomax[i]*norms[i] for i in range(len(norms))])
        return self.photoflux[band]
    
    def calc_continuum(self,wlength,fl,t):
        flux=copy(fl)
        photometry=self.get_photometry()
        
        BmV=np.power(mag_base,photometry['B-V'](t))
        Vmr=np.power(mag_base,photometry['V-r'](t))
        rmi=np.power(mag_base,photometry['r-i'](t))
        
        flux+=1.#add back the offset
        
        
        
        eq1=[IXn(wlength,flux,'V',2)-BmV*IXn(wlength,flux,'B',2),
             IXn(wlength,flux,'V',1)-BmV*IXn(wlength,flux,'B',1),
             IXn(wlength,flux,'V',0)-BmV*IXn(wlength,flux,'B',0)]
        
        eq2=[IXn(wlength,flux,'r',2)-Vmr*IXn(wlength,flux,'V',2),
             IXn(wlength,flux,'r',1)-Vmr*IXn(wlength,flux,'V',1),
             IXn(wlength,flux,'r',0)-Vmr*IXn(wlength,flux,'V',0)]
    
        eq3=[IXn(wlength,flux,'i',2)-rmi*IXn(wlength,flux,'r',2),
             IXn(wlength,flux,'i',1)-rmi*IXn(wlength,flux,'r',1),
             IXn(wlength,flux,'i',0)-rmi*IXn(wlength,flux,'r',0)]
        LHS=np.array([eq1,eq2,eq3])
        
        RHS=np.array([BmV*IXn(wlength,flux,'B',3)-IXn(wlength,flux,'V',3),
                      Vmr*IXn(wlength,flux,'V',3)-IXn(wlength,flux,'r',3),
                      rmi*IXn(wlength,flux,'r',3)-IXn(wlength,flux,'i',3)])
        b,c,d=np.linalg.solve(LHS,RHS)
        
        
        continuum=wlength**3+b*wlength**2+c*wlength+d
       
        
        return continuum
    
    def timeseries(self,offset=1.,lines=[10500.,11000.]):
        phase,wlength,flux=self.extract_data()
        plot_timeseries(phase, wlength, flux, offset, lines) 
        
    def multi_timeseries(self,fname):
        try:
            os.mkdir(fname)
        except:
            pass
        phase,wlength,flux=self.extract_data()
        limits=np.argmin(np.abs(phase))#get xlim and ylim from the closest to peak
        xmin,xmax=wlength[limits][0],wlength[limits][-1]
        xmin,xmax=7800,8800
        ymin,ymax=min(flux[limits]),max(flux[limits])
        dy=ymax-ymin
        ymin-=dy
        ymax+=dy
        for i in range(len(phase)):
            p.figure()
            p.plot(wlength[i],flux[i])
            p.ylim(ymin,ymax)
            p.xlim(xmin,xmax)
            p.title('t='+str(phase[i]))
            p.xlabel("wavelength (A)", fontsize=20)
            p.ylabel("flux + offset",fontsize=20)
            p.savefig(fname+'/%04d.png'%i)
            p.close()
        os.chdir(fname)
        os.system('ffmpeg -r %d -i %s.png -vcodec libx264 -y %s.mp4 -qscale 0' %(10,'%04d',fname))#create video
    
        os.chdir('..')#copy video to main directory and delete files

        
        
        
class Sav_Template(Template):
    def __init__(self,name):
        self.name='meanspec'+name#files start with 'meanspec'
        self.get_files() 
        self.continuum='one'
        self.extrapolate_lightcurve=False
        
        self.photometry_data='colormeans_corrected.pkl'
        self.light_curve_data='light_curves.dat'
        
    def get_files(self):
        out={}
        for fname in os.listdir(SAV_FOLDER):
            if fname.startswith(self.name) and fname.endswith('.sav'):
                try:
                    phase=float(fname[len(self.name)+1:-len('.sav')])#add 1 because expecting _
                    out[phase]=SAV_FOLDER+fname
                except ValueError:#the file has a longer name that contains the required name. e.g. Icbroad
                    pass
        self.files=out
        return self.files
    
    def extract_data(self):
        #extracts the data and returns it in a more usable format for plotting
        try:
            return [self.phase,self.wlength,self.flux] #check to see if it has already been run
        except AttributeError:
            phase_keys=np.sort(self.files.keys())
            phase=[]
            wlength=[]
            flux=[]
            error_p=[]
            error_n=[]
            maxx=[]
            minn=[]
            for i in phase_keys:
                sav=readsav(self.files[i])
                nonzero=np.nonzero(sav.fmean)[0]#remove zero padding
                try:
                    lo=nonzero[0]
                    hi=nonzero[-1]+1#need to add 1 due to how slices work
                except IndexError:#all fluxes are zero so don't use this phase
                    continue             
                
                phase.append(i)
                wlength.append(sav.wlog[lo:hi])
                flux.append(sav.fmean[lo:hi])
                maxx.append(sav.fmean[lo:hi]+sav.fmeanp[lo:hi])
                minn.append(sav.fmean[lo:hi]-sav.fmeann[lo:hi])
                error_p.append(sav.fsdev[lo:hi])
                error_n.append(sav.fsdev[lo:hi])
                
                
                #l,f=pload(self.files[i])
                #wlength.append(l)
                #flux.append(f)
            self.phase=np.array(phase)
            self.wlength=np.array(wlength)
            self.flux=np.array(flux)
            self.error_p=np.array(error_p)
            self.error_n=np.array(error_n)
            self.max_flux=np.array(maxx)
            self.min_flux=np.array(minn)
            return [self.phase,self.wlength,self.flux]
 
 
class Pickle_Template(Template):
    def __init__(self,fname,continuum='zero'):
        self.photometry_data='colormeans_corrected.pkl'
        self.light_curve_data='light_curves.dat'
        self.data=pload(fname)
        self.continuum=continuum
        self.extrapolate_lightcurve=False
    
    def extract_data(self):
        #extracts the data and returns it in a more usable format for plotting
        try:
            return [self.phase,self.wlength,self.flux] #check to see if it has already been run
        except AttributeError:
            phase_keys=np.sort(self.data.keys())
            phase=[]
            wlength=[]
            flux=[]
            error_p=[]
            error_n=[]
            for i in phase_keys:
                try:
                    l,f,var=self.data[i]
                except ValueError:
                    l,f,var,_=self.data[i]#4th entry is continuum      
                
                phase.append(i)
                wlength.append(l)
                flux.append(f)
                error_p.append(var)
                error_n.append(var)
            self.phase=np.array(phase)
            self.wlength=np.array(wlength)
            self.flux=np.array(flux)
            self.error_p=np.array(error_p)
            self.error_n=np.array(error_n)
            return [self.phase,self.wlength,self.flux]


data_folder='template_data/'
light_curve_folder='photometry/'
smoothing_scale=100.
class lnw_Template(Template):
    def __init__(self,fname,continuum='one',photometry='default'):
        self.fname=data_folder+fname+'.lnw'
        self.continuum=continuum
        f=open(self.fname,'r')
        lines=f.readlines()
        f.close()
        _,_,minl,maxl,continuum_length,self.sn_name,_,self.sn_type,_,_=lines[0].split()
        self.I=int(continuum_length)+2
        self.phase_info=lines[self.I].split()[0]
        self.extrapolate_lightcurve=False
        if photometry=='default':
            photometry=fname
        if photometry:
            self.get_photo_data(light_curve_folder+photometry+'.csv')
    
    def extract_data(self):
        #extracts the data and returns it in a more usable format for plotting
        try:
            return [self.phase,self.wlength,self.flux] #check to see if it has already been run
        except AttributeError:
            f=open(self.fname,'r')
            lines=f.readlines()
            f.close()
            phase=[]
            wlength=[]
            flux=[]
            for p in lines[self.I].split()[1:]:
                phase.append(float(p))
                flux.append([])
            for line in lines[self.I+1:]:
                vals=line.split()
                wlength.append(float(vals[0]))
                for i in range(len(phase)):
                    flux[i].append(float(vals[i+1]))
            '''remove zero padding'''
            tflux=[]
            tphase=[]
            twlength=[]
            for i in range(len(phase)):
                nonzero=np.nonzero(flux[i])[0]#remove zero padding
                try:
                    lo=nonzero[0]
                    hi=nonzero[-1]+1#need to add 1 due to how slices work
                except IndexError:#all fluxes are zero so don't use this phase
                    continue 
                tphase.append(phase[i])
                tflux.append(np.array(flux[i][lo:hi]))
                twlength.append(np.array(wlength[lo:hi]))
                
            self.phase=np.array(tphase)
            self.flux=np.array(tflux)
            self.wlength=np.array(twlength)     
            return [self.phase,self.wlength,self.flux]

        
class flm_Template(Template):
    def __init__(self,sn_name,continuum='zero',photometry='default',extrapolate_lightcurve=False):
        self.extrapolate_lightcurve=extrapolate_lightcurve
        self.sn_name=sn_name
        self.continuum=continuum
        self.folder=data_folder+sn_name+'/'
        self.fft_files={}
        self.error_files={}
        self.spectra_files={}
        f=open(self.folder+'Spec.JD.phases','r')
        lines=f.readlines()
        f.close()
        for line in lines[2:]:
            fname,_,phase=line.split()
            p=float(phase)
            self.spectra_files[p]=self.folder+fname
            self.fft_files[p]=self.folder+fname+'-fft-residual'
            self.error_files[p]=self.folder+'variance/'+fname[:-6]+'-sigma-z.flm'
        
        
        if photometry=='default':
            photometry=sn_name
        if photometry:
            self.get_photo_data(light_curve_folder+photometry+'.csv')
    
    def extract_data(self):
        #extracts the data and returns it in a more usable format for plotting
        try:
            return [self.phase,self.wlength,self.flux] #check to see if it has already been run
        except AttributeError:
            phases=np.sort(self.spectra_files.keys())
            wlength=[]
            flux=[]
            for p in phases:
                f=open(self.spectra_files[p],'r')
                lines=f.readlines()
                f.close()
                tl=[]
                tf=[]
                for line in lines:
                    a,b=line.split()
                    tl.append(float(a))
                    tf.append(float(b))
                wlength.append(np.array(tl))
                flux.append(np.array(tf))
            self.phase=np.array(phases)
            self.wlength=np.array(wlength)
            self.flux=np.array(flux)
            self.get_errors()
            return [self.phase,self.wlength,self.flux]
        
    def get_errors(self):
        try:
            return [self.error_p,self.error_n] #check to see if it has already been run
        except AttributeError:
            phases,wlength,flux=self.extract_data()
            
            errors=[]
            for i in range(len(phases)):
                phase=phases[i]
                try:
                    f=open(self.error_files[phase],'r')
                    lines=f.readlines()
                    f.close()
                    lambdas=[]
                    var=[]
                    for line in lines:
                        a,b=line.split()
                        lambdas.append(float(a))
                        var.append(float(b))
                except IOError:
                    print 'no error information for phase %g so resorting to Fourier analysis' %phase
                    f=open(self.fft_files[phase],'r')
                    lines=f.readlines()
                    f.close()
                    lambdas=[]
                    tvar=[]
                    for line in lines:
                        a,b,_=line.split()
                        lambdas.append(float(a))
                        tvar.append(float(b))
                    var=[]
                    for l in lambdas:
                        lo,hi=reduce_range(lambdas,l-smoothing_scale/2.,l+smoothing_scale/2.,indices=True)
                        var.append(np.std(tvar[lo:hi]))
                except KeyError:
                    print 'Error! no data exists for phase %g. Using 10 percent error bars' %phase
                    lambdas=wlength[i]
                    var=0.1*np.abs(flux[i]+1)   
                error_spec=interp1d(lambdas,var,bounds_error=False,fill_value=var[-1])    
                errors.append(error_spec(wlength[i]))
                
            self.error_p=np.array(errors)
            self.error_n=np.array(errors)
            return [self.error_p,self.error_n]
                


            

def IXn(wlength,flux,band,n):
    X=get_band(band)
    if min(wlength)>X.minl or max(wlength)<X.maxl:#data doesn't cover entire band so must normalise
        lmin=min(min(wlength),X.minl)
        lmax=max(max(wlength),X.maxl)
        denom=integrate(lambda z:X(z)*(z**n),lmin,lmax, N=100)
        a=integrate(lambda z:X(z)*(z**n),X.minl,lmin, N=100)
        b=integrate(lambda z:X(z)*(z**n),lmax,X.maxl, N=100)
        norm=1.+(a+b)/denom
    else:
        norm=1.
    y=flux*X(wlength)*(wlength**n)
    return norm*integrate_data(wlength,y)
 
def calc_photospectra(phase,wlength,flux,lc_phase,lc_flux,band='V'):
    lc=interp1d(lc_phase,lc_flux,bounds_error=False,fill_value=0)#this will take no power from spectra without lc data
    band_f=get_band(band)
    out=[]
    norms=[]
    for i in range(len(phase)):
        norm=integrate_data(wlength[i], flux[i]*band_f(wlength[i]))
        norm=lc(phase[i])/norm
        out.append(flux[i]*norm)
        norms.append(copy(norm))
    return (np.array(out),norms)               
    
    
def apply_window(source,window_data,ts,lambdas,band='V',minmax=False):#this will fail if you try to run photomax before mean
    window=interp1d(ts,window_data,bounds_error=False,fill_value=0.)
    phases,wlength,_=source.extract_data()
    if minmax=='min':
        flux=source.photomin[band]
    elif minmax=='max':
        flux=source.photomax[band]
    elif minmax=='extreme_max':
        flux=[(source.max_flux+1.)[i]*source.photonorms[band][i] for  i in range(len(source.photonorms[band]))]#assumes 1 for continuum
    elif minmax=='extreme_min':
        flux=[(source.min_flux+1.)[i]*source.photonorms[band][i] for  i in range(len(source.photonorms[band]))]
    else:
        flux=source.get_photospectra(band)
    interps=[interp1d(wlength[i],flux[i],bounds_error=False) for i in range(len(phases))]
    
    fluxes=[]
    for i in range(len(phases)):
        fluxes.append(interps[i](lambdas)*window(phases[i]))
        
    fluxes=np.array(fluxes).T
    out=[]
    for i in range(len(fluxes)):
        indices=np.isnan(fluxes[i])#filters out the nan values
        out.append(integrate_data(phases[~indices],fluxes[i][~indices]))
    return out    
        
def plot_timeseries(phase,wlength,flux,offset=0.05*1e-12,lines=[10500.,11000.]):
    colourdict=linear_gradient('#ff0000', '#3c26b7', len(phase))
    for i in range(len(phase)):
        p.plot(wlength[i],np.array(flux[i])+offset*i,colourdict['hex'][i])
        p.plot(lines,[offset*i]*2,colourdict['hex'][i])
        p.text(lines[0]+(2*(i%2)-1)*(lines[1]-lines[0]), offset*i, phase[i], fontsize=12)
    p.xlabel("wavelength (A)", fontsize=20)
    p.ylabel("flux + offset",fontsize=20)

def get_band(band):#returns the band as a function by reading the data
    fname='band_pass_filters/%s.txt'%band
    f=open(fname,'r')
    wlength=[]
    flux=[]
    for line in f:
        a,b=line.split()
        wlength.append(float(a))
        flux.append(float(b))
    f.close()
    #sort the list of fluxes
    together=zip(wlength,flux)
    sorted_together=sorted(together)
    
    wlength=np.array([x[0] for x in sorted_together])
    flux=np.array([x[1] for x in sorted_together])
    i=0
    while flux[i]<0.001:
        i+=1
    minl=wlength[i-1]
    while flux[i]>0.001:
        i+=1
    maxl=wlength[i]        
    out=Band(interp1d(wlength,flux,bounds_error=False,fill_value=0.),minl,maxl)#band template as a function. 0. outside support
    return out


def get_photo(data,colour):#returns the band as a function by reading the data
    ts,mags,errors=data
    indices=np.isnan(mags)
    out_t=np.array(ts)[~indices]
    out_mag=mags[~indices]
    ts=np.array(out_t)
    mags=np.array(out_mag)+zeropoints[colour]
    return lambda x:np.interp(x,copy(ts),copy(mags))

'''
def smooth_spectrum(lambdas,flux,var,vexp=0.00333,nsig=5,uniform=True):#shouldn't need nsig or uniform
    lambdas=np.array(lambdas)
    flux=np.array(flux)
    var=np.array(var)
    I=np.ones_like(lambdas)
    li=np.outer(lambdas,I)
    lj=np.outer(I,lambdas)
    epsj=np.outer(I,var)
    sigi=li*vexp
    fj=np.outer(I,flux)
    W=(1./epsj)*np.exp(-((li-lj)/sigi)**2)
    W=np.nan_to_num(W)#get rid of nans
    Wf=W*fj
    out=np.sum(Wf,0)/np.sum(W,0)
    return out
'''

'''
def smooth_spectrum(lambdas,flux,var,vexp=0.00333/2.,nsig=5,uniform=True):#this version is vectorised
    uniform=False
    if uniform:#assume a uniform spacing in lambda
        dl=lambdas[1]-lambdas[0]#same for all lambdas by assumption
        sigmas=vexp*lambdas
        nweights=int(nsig*max(sigmas)/dl)#number (+ and -) of wavelengths used to smooth
        padding=np.zeros(nweights)
        invar=1./var
        invar=np.concatenate((padding,invar,padding))
        flux=np.concatenate((padding,flux,padding))#pad with zeros        
        weights=[]
        weighted_flux=[]
        for i in range(-nweights,nweights+1):
            W=np.exp(-0.5*(i*dl/sigmas)**2)*invar[i+nweights:len(invar)+i-nweights]
            weights.append(copy(W))
            weighted_flux.append(copy(W*flux[i+nweights:len(invar)+i-nweights]))
        weights=np.array(weights)
        weighted_flux=np.array(weighted_flux)
        out=np.sum(weighted_flux,0)/np.sum(weights,0)
        return out
    else:
        out=[]
        for i in range(len(lambdas)):#not vectorised. try to vectorise if needs speeding up
            sigma=vexp*lambdas[i]
            gauss=gaussian(sigma,mu=lambdas[i])
            lo,hi=reduce_range(lambdas,lambdas[i]-nsig*sigma,lambdas[i]+nsig*sigma,indices=True)
            weights=gauss(lambdas[lo:hi])*(1./var[lo:hi])
            out.append(np.sum(weights*flux[lo:hi])/np.sum(weights))
        return np.array(out)
'''

def smooth_spectrum(lambdas,flux,var,vexp=0.00333/2.,nsig=5,uniform=True):
    out=[]
    for i in range(len(lambdas)):#not vectorised. try to vectorise if needs speeding up
        sigma=vexp*lambdas[i]
        gauss=gaussian(sigma,mu=lambdas[i])
        lambda_range=min(lambdas[i]-lambdas[0],lambdas[-1]-lambdas[0],nsig*sigma)
        lo,hi=reduce_range(lambdas,lambdas[i]-lambda_range,lambdas[i]+lambda_range,indices=True)
        if lo==hi:
            out.append(flux[lo])
        else:
            weights=gauss(lambdas[lo:hi])*(1./var[lo:hi])
            out.append(np.sum(weights*flux[lo:hi])/np.sum(weights))
    return np.array(out)

def get_min_cubic_spline(wlength,flux):
    '''finds the minimum of a feature but fitting a cubic spline to the nearest N points
    around the smallest value, where N is given by the width argument''' 
    
    
    spl,X=fit_spline(wlength,flux)
    
    if VERBOSE:
        x,y=plot_spline(spl,X)
        p.plot(x,y)
    
    y_min=np.inf
    for i in range(len(spl)):
        a,b,c,d=spl[i]
        x0=(-c+np.sqrt(c**2-3.*d*b))/(3.*d)+X[i]#minimum of the cubic 
        if x0<=X[i] and x0<=X[i+1]:#check that it is in the range of validity of the cubic
            y0=a+b*x0+c*x0**2+d*x0**3
            if y0<y_min:
                y_min=y0
                out=x0
    return (out,'l')#for backwards compatibility. Not a useful fit

def get_min_polyfit_quadratic(wlength,flux):
    '''finds the minimum of a feature but fitting a cubic function'''
    offset=wlength.mean()#removing this offset improves the quality of the fit
    scale=wlength.std()
    x=(wlength-offset)/scale
    
    c,b,a=np.polyfit(x,flux,2)
    if c>0:
        x0=-(scale*b)/(2*c)+offset
    else:#quadratic has no minimum, only max
        x0=np.nan
        
    test=(wlength[-1]-offset)/scale
    lr='l' if c*(test**2)+b*test>0 else 'r'
    if VERBOSE:
        x=np.linspace(wlength[0],wlength[-1],500)
        xx=(x-offset)/scale
        y=a+b*xx+c*(xx**2)
        p.plot(x,y)
        x1=(x0-offset)/scale
        p.plot(x0,a+b*x1+c*(x1**2),'bo')
    
        
    return (x0,lr)


def get_min_polyfit_cubic(wlength,flux):
    #finds the minimum of a feature but fitting a cubic function
    offset=wlength[0]#removing this offset improves the quality of the fit
    d,c,b,a=np.polyfit(wlength-offset,flux,3)
    x0=(-c+np.sqrt(c**2-3.*d*b))/(3.*d)+offset
    
    test=wlength[-1]-offset
    lr='l' if d*(test**3)+c*(test**2)+b*test>0 else 'r'
    
    if VERBOSE:
        x=np.linspace(wlength[0],wlength[-1],500)-offset
        y=a+b*x+c*(x**2)+d*(x**3)
        p.plot(x+offset,y)
        x1=x0-offset
        p.plot(x0,a+b*x1+c*(x1**2)+d*(x1**3),'bo')
    
            
    return (x0,lr)


def get_min_minimum(wlength,flux):
    i=np.nanargmin(flux)
    x0=wlength[i]
    if VERBOSE:
        p.plot(x0,flux[i],'bo')
    lr='l'#we will never be outside the allowed range so this should be irrelevent
    
            
    return (x0,lr)
    


def get_min_polyfit_quartic(wlength,flux,centre):
    '''finds the minimum of a feature but fitting a quartic function''' 
    offset=wlength[0]#removing this offset improves the quality of the fit
    e,d,c,b,a=np.polyfit(wlength-offset,flux,4)
    
    f=lambda x:b+2*c*x+3*d*(x**2)+4*e*(x**3)#dy/dx. want to equal zero
    df=lambda x: 2*c+6*d*x+12*e*(x**2)
    x0=newton(f,df,centre)+offset
    test=wlength[-1]-offset
    
    lr='l' if e*(test**4)+d*(test**3)+c*(test**2)+b*test>0 else 'r'
    
    if VERBOSE:
        x=np.linspace(wlength[0],wlength[-1],500)-offset
        y=a+b*x+c*(x**2)+d*(x**3)+e*(x**4)
        p.plot(x+offset,y)
        x1=x0-offset
        p.plot(x0,a+b*x1+c*(x1**2)+d*(x1**3)+e*(x1**4),'bo')
    
    
    return (x0,lr)

            
def get_min(wlength,mean,cov,nmc=False,fitting='cubic_poly',
            width=40,fixed_range=False,centre=None,full_range=True,smoothing=True,smoothing_matrix='none'):
    fitting_dict={'cubic_spline':get_min_cubic_spline,'cubic_poly':get_min_polyfit_cubic,'minimum':get_min_minimum,
                  'quartic_poly':lambda x,y:get_min_polyfit_quartic(x,y,centre),'quadratic_poly':get_min_polyfit_quadratic}
    global RECALC
    #cov should be either a covarience matrix or a vector of sigmas
    
    var1d=(len(cov.shape)==1)#true if 1D sigma supplied, false if 2D covarience
    
    find_min=fitting_dict[fitting]
    
    left,right=fixed_range
    lo,hi=reduce_range(wlength,left,right,include_ends=False,indices=True)
        
    l_min=copy(wlength[0])
    l_max=copy(wlength[-1])
    
    if var1d:
        sig=cov
    else:
        var=np.diagonal(cov)
        sig=np.sqrt(var) 
        try:
            L=np.linalg.cholesky(cov)#LL^T=cov
        except LinAlgError:#numerical issue, cannot cholesky decompose matrix
            print 'ERROR! covariance matrix is not positive definite. Resorting to diagonal components only.'
            var1d=True
           
    
    if not nmc:
        N=len(wlength[lo:hi])
        nmc=int(N*(np.log(N)**2))
    
    if VERBOSE:
        nmc=50#don't want to plot 3000 lines
        colourdict=linear_gradient('0000ff','ff0000',n=nmc)
        max_fails=1.1*nmc#when verbose want to plot every fit, even if they fail
    else:
        max_fails=0.25*nmc      
        
        
    diff=wlength[1:]-wlength[:-1]
    uniform=(max(diff)-min(diff))<np.mean(diff)*0.001
       
    failures=0
    MC_mins=[]
    l=lo
    r=hi
    for i in range(nmc):
        progress_bar(i,nmc)
        
        
        if var1d:
            tflux=np.random.normal(mean,sig)
        else:
            tflux=np.random.normal(np.zeros_like(mean),np.ones_like(mean))
            tflux=np.dot(L,tflux)+mean
        if smoothing:
            if smoothing_matrix=='none':
                tflux=smooth_spectrum(wlength,tflux,sig,uniform=uniform)#the l:r prevents measurement outside the desired range affecting the smoothing
                #tflux=np.concatenate(([0]*l,tflux,[0]*(len(wlength)-r)))#pad with zeros so it is the correct length
            else:
                tflux=mat_x_vec(smoothing_matrix,tflux)
        
        #if VERBOSE:
            #p.plot(wlength,tflux)
        if not full_range:
            minn=wlength[lo:hi][np.argmin(tflux[lo:hi])]
            l,r=reduce_range(wlength[lo:hi],minn-width,minn+width,indices=True)
            l+=lo
            r+=lo
        else:
            l,r=lo,hi
        x0,lr=find_min(wlength[l:r],tflux[l:r])
        if np.isnan(x0) or x0<l_min or x0>l_max: #remove for smallfit
            failures+=1
            if not full_range:
                x0=(wlength[l] if lr=='l' else wlength[r-1])
                MC_mins.append(x0)
        else:
            if x0<wlength[l] or x0>wlength[r-1]:
                if full_range:
                    failures+=1
                else:
                    x0=(wlength[l] if lr=='l' else wlength[r-1])#take minimum value within range of fit
            MC_mins.append(x0)
        if VERBOSE:
            #p.plot(wlength[l:r],usmooth[l:r]+i*5,color= colourdict['hex'][i],marker='x',ls='None')
            #p.plot(wlength[l:r],tflux[l:r]+i*5,color= colourdict['hex'][i])
            p.plot(wlength,tflux,color= colourdict['hex'][i])
            #p.plot(wlength,tflux,color= colourdict['hex'][i],marker='x',ls='None')
            p.axvline(x0)
                
        if failures>max_fails:
            print '\n'
            print 'The %s method does not work for this spectrum' %fitting
            return (np.nan,np.nan)
    print '\nthere were %d failures' %failures
    print '\n'
    MC_mins=np.array(MC_mins)
    
    return (np.median(MC_mins),np.percentile(MC_mins,15.865),np.percentile(MC_mins,84.135)) 
    
    
def get_regression(x,y,sig):
    '''uses regression to estimate the mean and the variance of the underlying distribution
    at the same positions as measured. Assumes gaussian kernel of width l'''
      
    #signal=10. #these were used when run in lambda space
    #l=100.
    
    signal=1.
    l=3.
    
    k=lambda a,b: signal*np.exp(-0.5*((a-b)/l)**2)#gaussian kernel
    
    var=sig**2
    
    x=np.array(x)
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
    
    ############################################################
    cov=0.5*(cov+cov.T)#cov should be symmetric, numerical factors might mess this up so restore it here
    return (np.array(mean).T[0],np.array(cov),K_xx*K_er)
    
    
def remove_continuum(x,y):#removes a cubic continuum
    a,b,c,d=np.polyfit(x,y,3)
    cont=a*(x**3)+b*(x**2)+c*x+d
    out=y/cont#divide out continuum
    return out
 
 
def rebin(x0,x1,y,cov,logsep=True,return_A=False):#re-bins a spectrum by averaging over N bins (todo, make general so that doesn't assume constant lambda difference)
    A=[]
    out_x=[]
    if logsep:
        x_edges=np.sqrt(x1[:-1]*x1[1:])
    else:
        x_edges=(x1[1:]+x1[:1])/2.
    x_edges=np.concatenate([[x0[0]],x_edges,[x0[-1]+1.]])#the +1 ensures that x0[-1] is still assigned to a bin
    for i in range(len(x1)):
        temp=np.ones_like(x0)*np.logical_and(x0>=x_edges[i],x0<x_edges[i+1])
        NN=np.sum(temp)
        if NN!=0:
            temp/=NN#normalise
            A.append(temp)
            out_x.append(x1[i])
    A=np.matrix(A)
    out_y=mat_x_vec(A,y)
    out_cov=A*cov*A.T#requires cov to be a matrix. todo: rewrite to allow 1d var
    out_x=np.array(out_x)
    if return_A:
        return (out_x,out_y,out_cov,A)
    else:
        return (out_x,out_y,out_cov)        
        
def rebin_set_bins(y,cov,N,return_A=False):#re-bins a spectrum by averaging over N bins (todo, make general so that doesn't assume constant lambda difference)
    Nbins=len(y)
    Nbins_new=int(np.ceil(float(Nbins)/N))
    A=[]
    for i in range(Nbins_new):
        temp=np.zeros(Nbins_new)
        temp[i]=1./N
        temp=np.repeat(temp,N)
        temp=temp[:Nbins]#this considers the case when the number of bins is not a multiple of N
        if i==Nbins_new-1:
            n=len(temp)-(Nbins_new-1)*N
            temp*=float(N)/float(n)#since the last bin is shorter, it needs a different normalisation
        A.append(temp)
    A=np.matrix(A)
    out_y=np.array(A*np.matrix(y).T).T[0]
    out_cov=A*cov*A.T#requires cov to be a matrix. todo: rewrite to allow 1d var
    if return_A:
        return (out_y,out_cov,A)
    else:
        return (out_y,out_cov)
        
    