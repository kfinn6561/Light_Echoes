'''
Created on Jul 22, 2015

@author: kieran
'''
from spectral_tools import lnw_Template, Pickle_Template
import pylab as p
import numpy as np
from colour_tools import kelly_colors
import os

fig=p.figure()
total=0
offset=0.1
colours =kelly_colors


sn=Pickle_Template('meanspec/meanspec_GPinterp.p',continuum='one')
sn.get_lightcurve()
phase, wlength,flux=sn.extract_data()#finds the range of phases and wavelthengths required
i=np.argmin(np.abs(phase))
p.plot(wlength[i],flux[i],color='k',lw=2,label='Template')


phase_info_dict={'0':'maximum light','1':'first observation','2':'unknown'}
for fname in os.listdir('template_data'):
    
    sn=lnw_Template('template_data/'+fname)
    
    if sn.sn_type=='IIb' and sn.phase_info=='0':
    #if sn.sn_name in ('sn1998fa','sn1993J','sn2003bg','sn2011ei'):
        print '%s is a type %s with data taken from %s.' %(sn.sn_name,sn.sn_type,phase_info_dict[sn.phase_info])
        total+=1
        phase,wlength,flux=sn.extract_data()
        i=np.argmin(np.abs(phase))
        p.plot(wlength[i],flux[i],color=colours[total-1],label=sn.sn_name)
        
        ''''p.plot(wlength[0],flux[0]+offset*phase[0],color=colours[total-1],label=sn.sn_name)
        for i in range(1,len(phase)):
            p.plot(wlength[i],flux[i]+offset*phase[i],color=colours[total-1])
        '''
print 'there are %d type IIb supernovae with phase info.' %total
p.legend()
p.show()







