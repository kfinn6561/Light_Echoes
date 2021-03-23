'''
Created on Sep 24, 2015

@author: kieran
'''
from general_tools import pload,pdump
import numpy as np

casa=pload('plot_data/casa_velocity_list.pkl')
individual=pload('plot_data/individual_velocity_list.pkl')

f=open('plots/velo_deluxetable.tex','w')

lines=['He','H_alpha','Ca_III']
les=['le3923','le2521','le2116']
sne=['Cas A','Template','sn2003bg','sn2011ei','sn1993J','sn1996cb','sn2008ax','sn2011dh']
le_plotnames={'le3923':'NE','le2521':'NW','le2116':'SW','unmodified':'w(t)=1'}

data={line:{le:{} for le in les} for line in lines}

for tdict in [casa,individual]:
    for line in tdict.keys():
        for le in tdict[line].keys():
            for sn in tdict[line][le].keys():
                data[line][le][sn]=tdict[line][le][sn]
pdump(data,'plot_data/velocity_list.pkl')
line_dict={'He':r'He I','H_alpha':r'H$\alpha$','Ca_III':'Ca II','Fe_II':'Fe 5169'}

f.write('\\tablehead{')
ll='&'
for sn in sne:
    ll+='&\\colhead{%s}' %sn
f.write(ll+'\n}\n')
f.write('\\startdata\n')
for line in lines:
    f.write(line_dict[line]+'&'+'&'*len(sne)+'\\\\\n')
    for le in les:
        ll='&%s'%le_plotnames[le]
        for sn in sne:
            try:
                x0,xmin,xmax=data[line][le][sn]
                xerr=(xmax-xmin)/2.
                if np.isnan(x0):
                    ll+='&'
                else:
                    ll+='&%.2f (%.2f)' %(x0,xerr)
            except KeyError:
                ll+='&'
        f.write(ll+'\\\\\n')
f.write('\\enddata\n')

f.close()
print 'done'