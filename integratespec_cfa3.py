#!/usr/bin/env python
import re,sys,string,math,os,types,exceptions,time,fcntl,shutil
# put the tools directory into the path
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS'])
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/SM')
sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/LESPEC')
import pylab as matlib
from   matplotlib.ticker import FormatStrFormatter
from texttable import txttableclass 
from sntable import sntableclass, intspecfilename, SPECTEMPLDIR, SPECREDUCTIONTYPES
from spectra import spectraclass 
import optparse


class intspecclass(txttableclass):
    def __init__(self,specreductiontype=None,tmax_MJD=None,z=None):
        txttableclass.__init__(self)
        # MLCS lc in restframe
        self.lc      = txttableclass()
        self.efflc   = txttableclass()
        self.intspec = txttableclass()
        self.rootdir = SPECTEMPLDIR
        self.snid=''

        self.lckey4peak = None
        self.peakMag    = None
        self.peaktime   = None

        # if self.refmag!=None: renormalize so that the peak mag is
        # self.refmag. The input spectra are normalized so that they
        # agree with the observed lc. Thus we need the observed peak
        # mags!
        # self.refmag = -20.0
        # MLCS peakmags in restframe
        # self.lcpeakmags_restframe = txttableclass()

        self.lambdamin  = self.lambdamax = None  

        # define the reference filter
        self.lcfilter = 'V'

        # define the reference filter and restphase column for the effective lightcurve.
        self.efflcfilter = 'V'
        self.efflcrestphase = 'restphase'

        # Define which type of spectra reduction to use (see SPECREDUCTIONTYPES in sntable for allowed types)
        self.specreductiontype = specreductiontype
        
        # what is the MJD at max?
        self.tmax_MJD = tmax_MJD
        self.z = z
        
    def loadspeclist(self,snid):
        listfilename = '%s/%s/%s.specinfo.mod' % (self.rootdir,snid,snid)
        if self.loadfile(listfilename):
            print 'ERROR! %s does not exist!!!!!' % listfilename
            sys.exit(0)
        self.configcols(['filename_used'],'s',visible=1)
        self.configcols(['restphase'],'f','%.2f',visible=1)
        self.configcols(['specweight','relweight','specweight_tot','relweight_tot'],'f','%f',visible=1)
        self.configcols(['w_lc2peak','wf'],'f','%f',visible=1)
        self.configcols(['minlambda','maxlambda'],'f','%.1f',visible=1)
        self.configcols(['minflux','maxflux'],'f','%.2e',visible=1)
        self.configcols(['mjdobs'],'f','%.2f')
        self.configcols(['spec','lambdakeys'],'o')
        self.configcols(['skip'],'d','%1d',visible=1)
        self.snid = snid

        speckeys = self.rowkeys()

        skiplistfilename =  '%s/%s/%s.skipspecnames' % (self.rootdir,snid,snid)    
        if os.path.isfile(skiplistfilename):
            skiplist = txttableclass()
            skiplist.loadfile(skiplistfilename)
            skipkeys = skiplist.rowkeys()
            self.speckeys = []
            for speckey in speckeys:
                self.setentry(speckey,'skip',0)
                for skipkey in skipkeys:
                    if self.getentry(speckey,'fname') == skiplist.getentry(skipkey,'filename'):
                        self.setentry(speckey,'skip',1)
                        break
                if self.getentry(speckey,'skip')==0:
                    self.speckeys.append(speckey)
        else:
            self.speckeys = speckeys
            self.setcol2value('skip',0)
        #self.printtxttable()
        #print skiplistfilename
        #sys.exit(0)
        self.speckeys = self.sortkeysbycols(self.speckeys,'restphase',asstring=0)
                
    def loadlc_MLCS(self,weightfunctionfilename=None,lefield=None,effectivelightcurvefilename=None,effmagcol=None,effrestphasecol=None):
        
        # load the lc file
        lcfilename = '%s/%s/%s.rest.interp.lc' % (self.rootdir,self.snid,self.snid)            
        if self.lc.loadfile(lcfilename):
            print 'ERROR! %s does not exist!!!!!' % lcfilename
            sys.exit(0)
        print 'light curve file %s loaded' % lcfilename

        self.lc.configcols(['mjdobs','restphase'],'f','%.4f')
        self.lc.configcols([self.lcfilter],'f','%.3f')
        self.lc.configcols(['speckeyminus','speckeyplus'],'d',visible=1)
        self.lc.configcols(['wminus_t','wminus_lc','wminus','wplus_t','wplus_lc','wplus'],'f','%.6f',visible=1)
        self.lc.configcols(['w_lc2peak','wf'],'f','%f',visible=1)

        # initialize the spline fit of the lightcurve
        self.lc.initspline('restphase',self.lcfilter,interpolationtype='linear')
        self.peakMag  = self.lc.spline(0.0)
        print 'peak Mag:',self.peakMag

        # get the keys of the spectra before and after a given epoch
        lckeys = self.lc.rowkeys()
        for lckey in lckeys:
            restphase = self.lc.asfloat(lckey,'restphase')
            minuskey = None
            pluskey = None
            for speckey in self.speckeys:
                #print restphase,self.getentry(speckey,'jd'),pluskey,minuskey
                if pluskey==None:
                    if self.asfloat(speckey,'restphase')>restphase:
                       pluskey = speckey
                if self.getentry(speckey,'restphase')<=restphase:
                    minuskey = speckey
            self.lc.setentry(lckey,'speckeyminus',minuskey)
            self.lc.setentry(lckey,'speckeyplus',pluskey)

            # calculate the ratio between peak and flux at epoch
            lcweight = math.pow(10,-0.4*(self.peakMag - self.lc.asfloat(lckey,self.lcfilter)))
            self.lc.setentry(lckey,'w_lc2peak',lcweight)

        if weightfunctionfilename!=None:
            if string.lower(weightfunctionfilename) == 'auto':
                if (lefield == None):
                    raise RuntimeError,"--weightfunctionfilename can only be 'auto' if --lefield is specified"
                tmp = re.sub('_.*','',lefield)
                weightfunctionfilename = '%s/%s/%s/%s.%s.subslitinfo.wf' % (options.wfdir,tmp,lefield,lefield,self.snid)
                print weightfunctionfilename
            if not os.path.exists(weightfunctionfilename):
                print 'ERROR! %s does not exist!!!!!' % weightfunctionfilename
                sys.exit(0)                
            if self.efflc.loadfile(weightfunctionfilename,takelastheader=1):
                print 'ERROR! %s does not exist!!!!!' % weightfunctionfilename
                sys.exit(0)
            print 'light curve weight function file %s loaded' % weightfunctionfilename
            if effrestphasecol!=None:
                self.efflcrestphase=effrestphasecol
            self.efflc.configcols([self.efflcrestphase,'dt','wfnorm1'],'f','%.4f')
            for wfkey in self.efflc.allrowkeys:
                self.efflc.setentry(wfkey,self.efflcrestphase,self.efflc.getentry(wfkey,'dt')*365.242199)
                                    
            # initialize the spline fit of the weight function
            (low_slope,low_slope_err,low_offset,low_offset_err)     =  self.efflc.straightline(self.efflc.allrowkeys[:4] ,self.efflcrestphase,'wfnorm1')
            (high_slope,high_slope_err,high_offset,high_offset_err) =  self.efflc.straightline(self.efflc.allrowkeys[-4:],self.efflcrestphase,'wfnorm1')
            print low_slope,low_slope_err,low_offset,low_offset_err
            print high_slope,high_slope_err,high_offset,high_offset_err
            
            self.efflc.initspline(self.efflcrestphase,'wfnorm1',interpolationtype='linear',low_slope=low_slope,high_slope=high_slope)
            for lckey in lckeys:
                wf = self.efflc.spline(self.lc.asfloat(lckey,'restphase'))
                if wf<0.0:
                    wf=0.0
                self.lc.setentry(lckey,'wf',wf)

        elif  effectivelightcurvefilename != None:
            if effmagcol!=None:
                self.efflcfilter=effmagcol
            if effrestphasecol!=None:
                self.efflcrestphase=effrestphasecol

            if self.efflc.loadfile(effectivelightcurvefilename):
                print 'ERROR! %s does not exist!!!!!' % effectivelightcurvefilename
                sys.exit(0)
            print 'effective light curve file %s loaded' % effectivelightcurvefilename
                
            if not self.efflcfilter in self.efflc.cols:
                print 'ERROR: the magcol of the effective lightcurve ',self.efflcfilter,' does not exist, only',self.efflc.cols
                sys.exit(0)
            if not self.efflcrestphase in self.efflc.cols:
                print 'ERROR: the restphase column of the effective lightcurve ',self.efflcrestphase,' does not exist, only',self.efflc.cols
                sys.exit(0)
                
            self.efflc.configcols([self.efflcrestphase],'f','%.4f')
            self.efflc.configcols([self.efflcfilter],'f','%.3f')

            self.efflc.initspline(self.efflcrestphase,self.efflcfilter,interpolationtype='linear')
            self.effpeakMag  = self.efflc.spline(0.0)
            print 'eff peak Mag:',self.effpeakMag
            
            for lckey in lckeys:
                efflcweight = math.pow(10,-0.4*(self.effpeakMag - self.efflc.spline(self.lc.asfloat(lckey,'restphase'))))
                self.lc.setentry(lckey,'wf', self.lc.getentry(lckey,'w_lc2peak')/efflcweight)
        else:
            for lckey in lckeys:
                self.lc.setentry(lckey,'wf',1.0)
            
        # get the peakmags
        #if self.refmag!= None:
        #    filename = '%s/%s/%s.MLCS.peakmags_restframe.dat' % (self.rootdir,self.snid,self.snid)            
        #    self.lcpeakmags_restframe.loadfile(filename)
        #    self.lcpeakmags_restframe.configcols(['U','B','V','R','I'],'f','%.3f',visible=1)
        #    if len(self.lcpeakmags_restframe.allrowkeys) != 1:
        #        raise RuntimeError,'Something is wrong with %s' % filename
            
        #self.lc.printtxttable(keys=self.lc.allrowkeys[:150])
        #sys.exit(0)
        return(0)

    
            
    def calcweights(self,phasemin=None,phasemax=None):
        lckeys = self.lc.rowkeys()
        self.setcol2value('specweight',0.0)
        self.setcol2value('specweight_tot',0.0)
        
#        self.lc.configcols(['w_lc2peak'],'f','%f',visible=1)

        for lckey in lckeys:
            

            restphase = self.lc.asfloat(lckey,'restphase')
            minuskey = self.lc.asint(lckey,'speckeyminus')
            pluskey = self.lc.asint(lckey,'speckeyplus')

            if phasemin!=None and restphase<phasemin:
                continue

            if phasemax!=None and restphase>phasemax:
                continue

            #lcweight = math.pow(10,-0.4*(self.peakMag - self.lc.asfloat(lckey,self.lcfilter)))
            #self.lc.setentry(lckey,'w_lc2peak',lcweight)

            # calculate the fracton of the two spectra contributing to that particular epoch
            if minuskey==None:
                restphaseplus = self.asfloat(pluskey,'restphase') 
                self.lc.setentry(lckey,'wminus_t',None)
                self.lc.setentry(lckey,'wminus_lc',None)
                self.lc.setentry(lckey,'wplus_t',1.0)
                self.lc.setentry(lckey,'wplus_lc',10**(-0.4*(self.lc.asfloat(lckey,self.lcfilter)-self.lc.spline(restphaseplus))))
                wminus = None
                wplus  = self.lc.getentry(lckey,'wplus_lc')*self.lc.getentry(lckey,'wf')
            elif pluskey==None:
                restphaseminus = self.asfloat(minuskey,'restphase') 
                self.lc.setentry(lckey,'wplus_t',None)
                self.lc.setentry(lckey,'wplus_lc',None)
                self.lc.setentry(lckey,'wminus_t',1.0)
                self.lc.setentry(lckey,'wminus_lc',10**(-0.4*(self.lc.asfloat(lckey,self.lcfilter)-self.lc.spline(restphaseminus))))
                wminus = self.lc.getentry(lckey,'wminus_lc')*self.lc.getentry(lckey,'wf')
                wplus  = None
            else:
                restphaseminus = self.asfloat(minuskey,'restphase') 
                restphaseplus = self.asfloat(pluskey,'restphase') 
            
                delta_restphase = restphaseplus-restphaseminus
                if delta_restphase<=0:
                    print 'ERROR!!!! delta_restphase=%f<=0! restphaseminus=%f  restphaseplus=%f' % (delta_restphase,restphaseminus,restphaseplus)
                    sys.exit(0)
                self.lc.setentry(lckey,'wminus_t',1.0-(restphase   - restphaseminus)/delta_restphase)
                self.lc.setentry(lckey,'wplus_t', 1.0-(restphaseplus - restphase   )/delta_restphase) 
                self.lc.setentry(lckey,'wminus_lc',10**(-0.4*(self.lc.asfloat(lckey,self.lcfilter)-self.lc.spline(restphaseminus))))
                self.lc.setentry(lckey,'wplus_lc',10**(-0.4*(self.lc.asfloat(lckey,self.lcfilter)-self.lc.spline(restphaseplus))))
                wminus = self.lc.getentry(lckey,'wminus_t')*self.lc.getentry(lckey,'wminus_lc')*self.lc.getentry(lckey,'wf')
                wplus  = self.lc.getentry(lckey,'wplus_t') *self.lc.getentry(lckey,'wplus_lc')*self.lc.getentry(lckey,'wf')
           
            self.lc.setentry(lckey,'wminus',wminus)
            self.lc.setentry(lckey,'wplus',wplus)
               
            if wminus != None:
                self.setentry(minuskey,'specweight',self.getentry(minuskey,'specweight')+wminus)
                #self.setentry(minuskey,'specweight_tot',self.getentry(minuskey,'specweight_tot')+wminus/self.lc.asfloat(lckey,'w_lc2peak'))
            if wplus != None:
                self.setentry(pluskey,'specweight',self.getentry(pluskey,'specweight')+wplus)
                #self.setentry(pluskey,'specweight_tot',self.getentry(pluskey,'specweight_tot')+wplus/self.lc.asfloat(lckey,'w_lc2peak'))

        #speckeys = self.rowkeys()
        # update the main table with phase, and normalization factors
        for speckey in self.speckeys:
            M = self.lc.spline(self.getentry(speckey,'restphase'))
            if self.specreductiontype in ['dered-warp','dered-galsub']:
                w_lc2peak =  math.pow(10,-0.4*(self.peakMag-M))
            elif self.specreductiontype == 'unmodified':
                raise RuntimeError,'"unmodified" needs to be checked if treatment of effective lightcurve works!!!!'
                w_lc2peak =  1.0
            elif self.specreductiontype == 'old':
                raise RuntimeError,'"old" is not supported with integratespec.py, use integratespec_old.py instead!'
            else:
                print 'ERROR: wrong specreductiontype=%s! only the following are allowed:',SPECREDUCTIONTYPES
                sys.exit(0)
            self.setentry(speckey,'w_lc2peak',w_lc2peak)       
            if  options.weightfunctionfilename!=None:
                wf = self.efflc.spline(self.getentry(speckey,'restphase'))
                if wf<0.0:
                    wf=0.0
                self.setentry(speckey,'wf',wf) 
            elif  options.effectivelightcurve!=None:
                effM = self.efflc.spline(self.getentry(speckey,'restphase'))
                w_efflc2peak =  math.pow(10,-0.4*(self.effpeakMag-effM))
                self.setentry(speckey,'wf',w_lc2peak/w_efflc2peak)  
            else:
                self.setentry(speckey,'wf',1.0)  

                
        fout=open('newlc.dat','w')
        self.lc.printtxttable(cols=['mjdobs','restphase',self.lcfilter,'speckeyminus','speckeyplus','wminus_t','wminus_lc','wminus','wplus_t','wplus_lc','wplus','w_lc2peak','wf'], file=fout)
        self.printtxttable()
        
    def loadspectra(self):
        speckeys = self.rowkeys()
        for speckey in speckeys:
            self.setentry(speckey,'spec',spectraclass())
            spec = self.getentry(speckey,'spec')
            if self.specreductiontype in ['dered-warp','dered-galsub']:
                filename = '%s/%s/%s-%s.norm' % (self.rootdir,self.snid,self.specreductiontype,self.getentry(speckey,'fname'))
            elif self.specreductiontype == 'unmodified':
                filename = '%s/%s/%s.norm' % (self.rootdir,self.snid,self.getentry(speckey,'fname'))
            elif self.specreductiontype == 'old':
                raise RuntimeError,'"old" is not supported with integratespec.py, use integratespec_old.py instead!'
            else:
                print 'ERROR: wrong specreductiontype=%s! only the following are allowed:',SPECREDUCTIONTYPES
                sys.exit(0)
            self.setentry(speckey,'filename_used',os.path.basename(filename))    
        
            if spec.loadfile(filename):
                print 'ERROR! %s does not exist!!!!!' % filename
                sys.exit(0)
               
            spec.configcols(['lambda','flux'],'f','%f')
            minlambda = spec.minentry('lambda')
            self.setentry(speckey,'minlambda',minlambda)
            maxlambda = spec.maxentry('lambda')
            self.setentry(speckey,'maxlambda',maxlambda)
            self.setentry(speckey,'minflux',spec.minentry('flux'))
            self.setentry(speckey,'maxflux',spec.maxentry('flux'))
        if self.lambdamin ==None:
            self.lambdamin = self.maxentry('minlambda',keys=self.speckeys)
        if self.lambdamax ==None:
            self.lambdamax = self.minentry('maxlambda',keys=self.speckeys)

    def getlambdakeys(self,speckey):
        lambdakeys = self.getentry(speckey,'spec').CUT_inrange('lambda',self.lambdamin,self.lambdamax)
        return lambdakeys
        
    def calcintegratedspectra(self):
        if (os.path.isfile('%s/%s/%s.lambdarange' % (self.rootdir,self.snid,self.snid))):
            lines=open('%s/%s/%s.lambdarange' % (self.rootdir,self.snid,self.snid)).readlines()
            self.lambdamin = float(lines[0])
            self.lambdamax = float(lines[1])

        # Determine if all the spectra have the same wavelength vector
        lambdadifferent_flag = False
        #speckeys = self.rowkeys()
        reflambdavals = self.getentry(self.speckeys[0],'spec').col_asfloat_list('lambda',keys=self.getlambdakeys(self.speckeys[0]))
        for speckey in self.speckeys:
            lambdakeys = self.getlambdakeys(speckey)
            self.setentry(speckey,'lambdakeys',lambdakeys)
            lambdavals = self.getentry(speckey,'spec').col_asfloat_list('lambda',keys=lambdakeys)
            if reflambdavals != lambdavals:
                print 'lambda vector of %s is different than reference %s' % (self.getentry(speckey,'fname'),self.getentry(self.speckeys[0],'fname'))
                lambdadifferent_flag = True

        self.intspec.configcols(['lambda'],'f',visible=1)
        self.intspec.configcols(['integrated_flux'],'f','%e',visible=1)

        # add up all the weights, then normalize so that sum of weights is 1.0, thus the V magnitude of the integrated spectra should be the same than the V band magnitude of the input spectra
        # NOTE: use column specweight for that since these are the weights with respect to the spectra that are normalized to the same magnitude
        sumweights = 0.0
        sumweights_tot = 0.0        
        for speckey in self.speckeys:
            sumweights     += self.getentry(speckey,'specweight')
            self.setentry(speckey,'specweight_tot',self.getentry(speckey,'specweight')/self.asfloat(speckey,'w_lc2peak'))
            sumweights_tot += self.getentry(speckey,'specweight_tot')        
        for speckey in self.speckeys:
            self.setentry(speckey,'relweight',self.getentry(speckey,'specweight')/sumweights)
            self.setentry(speckey,'relweight_tot',self.getentry(speckey,'specweight_tot')/sumweights_tot)
            
        # if the lambda vectors are not different: easy!!!! just add up the spectra!
        if not lambdadifferent_flag:
            print '####### Wavelength vectors are the same'
            print '####### Integrating spectra'
            irange = range(len(reflambdavals))
            for i in irange:
                sum = 0.0
                refspec      = self.getentry(self.speckeys[0],'spec')
                reflambdakey = self.getentry(self.speckeys[0],'lambdakeys')[i]
                reflambda    = refspec.getentry(reflambdakey,'lambda')
                for speckey in self.speckeys:
                    lambdakey = self.getentry(speckey,'lambdakeys')[i]
                    spec = self.getentry(speckey,'spec')
                    weightedflux = spec.getentry(lambdakey,'flux')*self.getentry(speckey,'specweight')
                    
                    # renormalize to a reference mag if wanted
                    #if self.refmag != None:                        
                    #    peakMag_restframe = self.lcpeakmags_restframe.getentry(self.lcpeakmags_restframe.allrowkeys[0],self.lcfilter)
                    #    weightedflux *= math.pow(10,-0.4*(self.refmag-peakMag_restframe))

                    sum +=  weightedflux
                    #print 'DDDD',sum,self.getentry(speckey,'spec').getentry(lambdakey,'flux'),self.getentry(speckey,'specweight')
                    if (reflambda != spec.getentry(lambdakey,'lambda')):
                        print 'ERROR: lambdas are not the same!!!!',reflambda,spec.getentry(lambdakey,'lambda')
                        sys.exit(0)
                #sum /= sumweights
                self.intspec.newrow({'lambda':  reflambda,'integrated_flux': sum})
        else:
            print '####### SPECTRA HAVE DIFFERENT WAVELENGTH VECTORS! initializing vectors so that they can be interpolated.'
            for speckey in self.speckeys:
                spec = self.getentry(speckey,'spec')
                spec.initializespectra(lambdacol='lambda',fluxcol='flux',
                                       accessmethod=None,lambdastep=2.0,
                                       interpolmethod='average',lambdaboxsize=2.0,verbose=1)
            lambdarange = range(int(self.lambdamin+1),int(self.lambdamax))
            print '####### Integrating spectra (interpolation)'
            
            for lambdaval in lambdarange:
                sum = 0.0
                for speckey in self.speckeys:
                    spec = self.getentry(speckey,'spec')
                    (errorflag,interpoltemplflux) = spec.interpolflux(lambdaval)
                    if errorflag:
                        print 'Error calculating flux for lambda=%f of file %s' % (lambdaval,self.getentry(speckey,'fname'))
                        sys.exit(0)
                    weightedflux = interpoltemplflux*self.getentry(speckey,'specweight')

                    sum += weightedflux
                #sum /= sumweights
                self.intspec.newrow({'lambda': lambdaval ,'integrated_flux': sum})
        
        if options.normalize2mag:
            lambdamin=5000
            lambdamax=6500
            sum = 0.0
            keys = self.intspec.allrowkeys
            for i in xrange(1,len(keys)-1):
                lambdaval = self.intspec.getentry(keys[i],'lambda')
                if lambdaval<lambdamin or lambdaval>lambdamax:
                    continue
                sum+=self.intspec.getentry(keys[i],'integrated_flux')*(self.intspec.getentry(keys[i+1],'lambda')-self.intspec.getentry(keys[i-1],'lambda'))*0.5
           
            print 'normalizing so that magnitude of flux in %f to %f is %f' % (lambdamin,lambdamax,options.mag4normalize)
            M = -2.5 * math.log10(sum)
            c = math.pow(10,-0.4*(options.mag4normalize - M)) 
            print 'normalizing factor: %.4e' % (c)
            for key in self.intspec.allrowkeys:
                self.intspec.setentry(key,'integrated_flux',self.intspec.getentry(key,'integrated_flux')*c)

    def saveintegratedspectra(self,outfilename=None,intspecdir=None,lefield=None):
        if outfilename == None:
            outfilename = intspecfilename(self.snid,self.specreductiontype,intspecdir=intspecdir,lefield=lefield)
            print outfilename
        outfilename = os.path.expanduser(outfilename)
        intspecdir = os.path.abspath(os.path.dirname(outfilename))
        if not os.path.isdir(intspecdir):
            os.makedirs(intspecdir)
        print 'Saving integrated spectrum into %s, wavelength range %f to %f' % (outfilename,self.lambdamin,self.lambdamax)
        self.intspec.save2file(outfilename)
        self.saveresults2file('%s.results' % outfilename)

    def saveresults2file(self,outfilename = None):
        if outfilename == None:
            if self.specreductiontype in SPECREDUCTIONTYPES:
                outfilename = '%s/%s/%s.%s.results' % (self.rootdir,snid,snid,self.specreductiontype)
        print 'Saving %s' % outfilename
        self.save2file(outfilename)
        
    def results2sntable(self,sntable,key):
        sntable.setentry(key,'Nspec',len(self.speckeys))
        sntable.setentry(key,'Nspecused',len(self.speckeys)-len(self.CUT_inrange('relweight',None,0.0)))
        
        sntable.setentry(key,'minlambda',self.lambdamin)
        sntable.setentry(key,'maxlambda',self.lambdamax)
        sntable.setentry(key,'relweightmax',self.maxentry('relweight_tot',keys=self.speckeys))
        sntable.setentry(key,'phasemin',self.minentry('restphase',keys=self.speckeys))
        sntable.setentry(key,'phasemax',self.maxentry('restphase',keys=self.speckeys))
        if sntable.getentry(key,'relweightmax')<0.25:
            sntable.setentry(key,'ok',1)
        elif sntable.getentry(key,'relweightmax')<0.35:
            sntable.setentry(key,'ok',2)
        elif sntable.getentry(key,'relweightmax')<0.45:
            sntable.setentry(key,'ok',3)
        else:
            sntable.setentry(key,'ok',0)
            
            

                        
if __name__=='__main__':

    parser = optparse.OptionParser(usage = 'integratespec.py sn1 sn2 sn3 ... [options]. This uses the lightcurve normalized spectra!!!!')
    parser.add_option('-t','--specreductiontype'  ,  default="dered-warp", type="string",
                      help='Define which type of spectra reduction to use (dered-warp, dered-galsub, unmodified, old)')
    parser.add_option('-s','--save'  , default=False, action="store_true",
                      help='save the integrated spectra and the results')
    parser.add_option('-n','--normalize2mag'  , default=False, action="store_true",
                      help='normalize the integrated spectra so that the convolution with R gives the instrumental magnitude specified with --mag4normalize')
    parser.add_option('-m','--mag4normalize'  ,  default=-8.0, type="float", 
                      help='if --normalize2mag, then normalize the integrated spectra so that the convolution with R gives the specified instrumental magnitude (default:-8)')
    parser.add_option('-i','--phasemin'  ,  default=None, type="float", 
                      help='specify the restphase minimum from which to integrate')
    parser.add_option('-a','--phasemax'  ,  default=None, type="float", 
                      help='specify the restphase maximum to which to integrate')
    parser.add_option('-c','--clearall'  ,   default=False, action="store_true",
                      help='Clears all fit results in sntable')
    parser.add_option('-o','--outfilename'  ,   default=None, type="string",
                      help='specify the output filename for the integrated spectra. (note: overwrites --outdir)')
    parser.add_option('--intspecdir'  , default=os.environ['LESPECINT_DIR'], type="string",
                      help='dir where the integrated spectra are saved (default=%default)')
    parser.add_option('--wfdir'  , default=os.environ['LESPECWF_DIR'], type="string",
                      help='dir where the window functions are (default=%default)')
    parser.add_option('-f','--lefield'  ,   default=None, type="string",
                      help='specify the light echo field, which is then put into the output filename')
    parser.add_option('-w','--weightfunctionfilename',   default=None, type="string",
                      help='specify the lightcurve weightfunction file. (\'auto\' possible)')
    parser.add_option('-l','--effectivelightcurve'  ,   default=None, type="string",
                      help='specify the effective lightcurve file.')
    parser.add_option('--effmagcol'  ,   default=None, type="string",
                      help='specify the magnitude column in the effective lightcurve file')
    parser.add_option('--effrestphasecol'  ,   default=None, type="string",
                      help='specify the restphase column in the effective lightcurve file')
    options, args = parser.parse_args()
    if not (options.specreductiontype in SPECREDUCTIONTYPES):
        print 'ERROR: wrong specreductiontype=%s! only the following are allowed:',SPECREDUCTIONTYPES
        sys.exit(0)

    if len(args)>0:
        if args[0].lower() == 'all':
            snlist = []
        else:
            snlist = args
    else:
        print 'List the SNe as arguments (or \'all\')'
        sys.exit(0)

    sntable = sntableclass()
    sntable.loadsntable(specreductiontype=options.specreductiontype)

    if options.clearall:
        sntable.setcols2value(['Nspec','Nspecused','relweightmax','minlambda','maxlambda','phasemin','phasemax','ok'],None)                    

    keys = sntable.getselectedkeys(snlist)    
    for key in keys:
        snid = sntable.getentry(key,'snid')
        print snid
        intspec = intspecclass(specreductiontype = options.specreductiontype,
                               tmax_MJD = sntable.getentry(key,'tmax_MJD'),
                               z = sntable.getentry(key,'z'))

        # load the list of the spectra. this list will be in self
        intspec.loadspeclist(snid)

        # load the lightcurve. this list will be in self.lc        
        if intspec.loadlc_MLCS(weightfunctionfilename=options.weightfunctionfilename,lefield=options.lefield,effectivelightcurvefilename=options.effectivelightcurve,effmagcol=options.effmagcol,effrestphasecol=options.effrestphasecol):
            sntable.setentry(key,'ok',-1)

            
            print 'WARNING: SKIPPING %s' % snid
            continue
        print 
        # calculate the weights for each spectra, based on the light-curve
        intspec.calcweights(phasemin=options.phasemin,phasemax=options.phasemax)
        intspec.loadspectra()
        
        intspec.calcintegratedspectra()
        if options.save:
            intspec.saveintegratedspectra(outfilename=options.outfilename,intspecdir=options.intspecdir,lefield=options.lefield)
            #intspec.saveresults2file()
        intspec.results2sntable(sntable,key)

        #intspec.lc.printtxttable()
#        intspec.printtxttable()
    #sntable.printtxttable()
    if options.save:
        sntable.savesntable()
        #sntable.setcol2value('relweightmax',0.0)
        #sntable.setcol2value('Nspec',0)
        #sntable.setcol2value('minlambda',0.0)
        #sntable.setcol2value('maxlambda',0.0)
        #sntable.setcol2value('phasemin',0.0)
        #sntable.setcol2value('phasemax',0.0)
        #sntable.setcol2value('ok',0)
        #sntable.save2file(sntable.filename,cols=['snid','source','typedesc','type','dm15','dm15_err','Nspec','relweightmax','minlambda','maxlambda','phasemin','phasemax','ok'])
