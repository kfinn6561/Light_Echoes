'''
Created on Jul 4, 2015

@author: kieran
'''

from general_tools import pload
import pylab as p
import numpy as np

lambdas=np.linspace(4800,9000,500)#this is set up for a specific run. Change if more general needed

def scale(x):#scale such that the average value is 1
    x=np.array(x)
    norm=np.sum(x)/len(x)
    return x/norm

def plot_compare(small,large,title):
    p.figure()
    p.plot(lambdas,scale(small))
    p.plot(lambdas,scale(large))
    p.legend(['Narrow Window','Full Window'])
    p.xlabel('Wavelength (A)')
    p.ylabel('Flux')
    p.title(title)
    p.savefig('comparison_%s.png' %title)
    

small_unmod,small_spectra=pload('small_window.dat')
large_unmod,large_spectra=pload('full_window.dat')


plot_compare(small_unmod,large_unmod,'Unmodified')
plot_compare(small_spectra[0], large_spectra[0], 'le2116')
plot_compare(small_spectra[1], large_spectra[2], 'le2521')
plot_compare(small_spectra[2], large_spectra[2], 'le3923')

p.show()