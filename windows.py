'''
Created on May 25, 2015

RUN SOURCE EXPORT.SOURCEME BEFORE STARTING


@author: kieran
'''
from spectral_tools import *
import pylab as p

supernova='sn1993j'
reduction_type='dered-warp'    

def top_hat(width,centre=0.5,norm=1.):#centre is fraction 0-1
    height=norm/width
    support=width*(np.array([0,1])-centre)
    return window_function(lambda x:height,support)

def triangle(left,right=False,norm=1.):
    if not right:
        right=left
    h=norm*2/(left+right)
    return window_function(lambda x:np.piecewise(x,[x<0,x>=0],[lambda s:h*(1+s/left),lambda s:h*(1-s/right)]),[-left,right])


print('Loading the data')
sp=spec_plotter()
sp.loadspeclist(supernova)
sp.specreductiontype=reduction_type
sp.loadspectra()

print 'plotting original timeseries'
p.figure()
sp.timeseries()
p.title('Original')
'''
print '\n'
print 'applying top hat window function'
x,y,z=sp.apply_window(top_hat(20))
p.figure()
plot_timeseries(x,y,z)
p.title('Top Hat')

print '\n\n'
print 'applying triangular window function'
x,y,z=sp.apply_window(triangle(10))
p.figure()
plot_timeseries(x,y,z)
p.title('Triangle')
'''
p.show()