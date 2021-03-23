'''
Created on Jul 22, 2015

@author: kieran
'''
from spectral_tools import *
from general_tools import pload,pdump
  
print('Loading the data')
sn1993j=spec_plotter()
sn1993j.loadspeclist('sn1993j')
sn1993j.specreductiontype='dered-warp'
sn1993j.loadspectra()
sn1993j.loadlc_MLCS()

sn=Template('IIb_1specperSN')

phase,flux=sn1993j.get_lightcurve()
colours=sn.get_photometry()

out={'V':flux}
out['B']=flux+colours['B-V'](phase)
out['r_sdss']=flux-colours['V-r'](phase)
out['i_sdss']=out['r_sdss']-colours['r-i'](phase)

pdump([phase,out],'light_curves.dat')
print 'done'