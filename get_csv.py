import numpy as np
import sys
from general_tools import pload,pdump

sn_name=sys.argv[1]


csv_name='photometry/%s.csv' %sn_name
pkl_name='photometry/%s.pkl' %sn_name

try:
    phase,flux=pload(pkl_name)
except IOError:
    new_name='photometry/%s_photo.pkl' %sn_name
    phase,flux=pload(new_name)
    pdump([phase,flux],pkl_name)
bands=flux.keys()

f=open(csv_name,'w')

st=''
for band in bands:
    st+=','+band
f.write(st)
for i in range(len(phase)):
    f.write('\n')
    st=str(phase[i])
    for band in bands:
        st+=','
        if not np.isnan(flux[band][i]):
            st+=str(flux[band][i])
    f.write(st)
f.close()
    