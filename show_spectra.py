'''
Created on Jul 21, 2015

@author: kieran
'''
from general_tools import pload
import pylab as p
from spectral_tools import *

print('Loading the data')
sn1993j=spec_plotter()
sn1993j.loadspeclist('sn1993j')
sn1993j.specreductiontype='dered-warp'
sn1993j.loadspectra()
sn1993j.loadlc_MLCS()

sn=Template('IIb_1specperSN')
sn.get_lightcurve(sn1993j)#this is obviously rubbish. Just a placeholder until I get the real data
phases,wlength,flux=sn.extract_data()#finds the range of phases required
photoflux=sn.get_photospectra()
_,wl,continua=pload('continua.pkl')

for i in range(len(phases)):
    p.figure()
    p.plot(wlength[i],photoflux[i])
    x=photoflux[i]/(continua[i]*(flux[i]+1))
    norm=(min(x)+max(x))/2
    p.plot(wl[i],continua[i]*norm)
    p.title(phases[i])
    
p.show()