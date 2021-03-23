import matplotlib as mpl
import pylab as pl
from pylab import rc
rc('axes', linewidth=1.2)
#from matplotlib.font_manager import FontProperties##


#FontProperties.set_weight('normal')
mpl.rcParams['font.size'] = 40.
#mpl.rcParams['font.size'] = 22.
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']#Computer Modern Roman']
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['text.usetex'] = True
#print mpl.rcParams['font.serif']
#mpl.rcParams['font.serif'] = 'Times New Roman'#Bitstream Vera Serif'
#print mpl.rcParams['font.serif']
mpl.rcParams['axes.labelsize'] = 22
mpl.rcParams['xtick.labelsize'] = 22.
mpl.rcParams['ytick.labelsize'] = 22.
#mpl.rocParams['axes.labelsize'] = 40
#mpl.rcParams['xtick.labelsize'] = 30.
#mpl.rcParams['ytick.labelsize'] = 30.
mpl.rcParams['xtick.major.size']= 10.
mpl.rcParams['xtick.minor.size']= 5.
mpl.rcParams['ytick.major.size']= 10.
mpl.rcParams['ytick.minor.size']= 5.

params = {'legend.fontsize': 40, #'legend.linewidth': 1, # legend.linewidth is not a valid rc parameter after I update my python packages for the software carpentry workshop at 225th AAS
          'legend.numpoints':1,
          'legend.handletextpad':1
      }
# check a list of valid parameters
# print pl.rcParams.keys()

pl.rcParams.update(params) 
