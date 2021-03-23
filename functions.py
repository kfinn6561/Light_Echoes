'''
Created on 14 Jan 2013

@author: Kieran Finn
'''
import numpy as np
#from pycuba import Vegas, Cuhre
from general_tools import progress_bar
from copy import copy

def Linspace(a,b,N):
    #a rewriting of numpy's linspace function that works for multi-dimensional arrays
    try:
        return np.linspace(a,b,N)#may want to add kwargs if becomes important
    except:
        c=np.linspace(0,1,N)
        e=np.ones_like(a)
        d=np.ones_like(c)
        return np.multiply.outer(d,a)+np.multiply.outer(d,b-a)*np.multiply.outer(c,e)

def simpson(f,a,b):
    #integrates the function f from a to b using simpson's rule
    return ((b-a)/6.)*(f(a)+4*f((a+b)/2)+f(b))


'''
integration_methods={'vegas':Vegas,'cuhre':Cuhre}
def Integrate(f,a,b,method='cuhre'):
    ndim=len(a)#a,b are ND vectors that contain the limits of the integration
    a=np.array(a)
    b=np.array(b)
    def Integrand(ndim, xx, ncomp, ff, userdata):
        # access the current parameters
        x=np.array([xx[i] for i in range(ndim.contents.value)])
        # compute the result
        x=a+(b-a)*x#because integral is over 0-1
        result = f(*x)
        # store the result (here only one component)
        ff[0] = result
        return 0
    out=integration_methods[method](Integrand,ndim)['results'][0]['integral']#this should be fine, but may need to generalise at some point
    out*=np.prod(b-a)#multiply by the jacobian
    return out
'''    
    
    

def integrate(f,a,b,N=10,method=simpson):
    x=Linspace(a,b,N)
    out=method(f,x[:-1],x[1:])
    return np.sum(out,0)

def integrate_nv(f,a,b,N=10,method=simpson,verbose=False):
    #same as integrate but not vectorised so slower but more stable
    x=Linspace(a,b,N)
    out=0.
    for i in range(len(x)-1):
        if verbose:
            progress_bar(i,N-1)
        out+=method(f,x[i],x[i+1])
    return out

def integrate_2d(f,x1,x2,y1,y2,Nx=10,Ny=10):
    #2 dimensional integral, written in vectorised form to speed things up.
    x=np.linspace(x1,x2,2*Nx-1)
    y=np.linspace(y1,y2,2*Ny-1)
    xl=x[:-2:2]
    xhf=x[1:-1:2]
    xh=x[2::2]
    yl=y[:-2:2]
    yhf=y[1:-1:2]
    yh=y[2::2]
    xl,yl=np.meshgrid(xl,yl)
    xhf,yhf=np.meshgrid(xhf,yhf)
    xh,yh=np.meshgrid(xh,yh)
    pre=(yh-yl)*(xh-xl)/36.
    out=pre*(f(xl,yl)+f(xl,yh)+4*f(xl,yhf)+f(xh,yl)+f(xh,yh)+4*f(xh,yhf)+4*f(xhf,yl)+4*f(xhf,yh)+16*f(xhf,yhf))
    return np.sum(out)
        

def integrate_data(x,y):
    #integrates a function, given a list of x and y values using trapezium rule
    x=np.array(x)
    y=np.array(y)
    h=x[1:]-x[:-1]
    out=0.5*(y[:-1]+y[1:])*h
    return np.sum(out)
    

def Heaviside(x): #Heaviside step function in a form that can be calculated for numpy arrays
    if x>0:
        return 1.
    else:
        return 0.

def HS(x):
    return np.piecewise(x, [x<0,x>=0], [0.,1.])

def top_hat(x,left,right):
    if type(x) is np.ndarray:
        return HS(x-left)-HS(x-right)
    else:
        return Heaviside(x-left)-Heaviside(x-right)
    
def log_n(x,base):#log to base n using numpy so will work on arrays
    '''log_n(x)=log_e(x)/log_e(n)'''
    return np.log(x)/np.log(base)


def fit_spline(X,Y):#fits a cubic spline to the data and returns the coefficients of the fit
    np1=len(X)
    n=np1-1
    a = Y[:]
    b = [0.0]*(n)
    d = [0.0]*(n)
    h = [X[i+1]-X[i] for i in xrange(n)]
    alpha = [0.0]*n
    for i in xrange(1,n):
        alpha[i] = 3/h[i]*(a[i+1]-a[i]) - 3/h[i-1]*(a[i]-a[i-1])
    c = [0.0]*np1
    L = [0.0]*np1
    u = [0.0]*np1
    z = [0.0]*np1
    L[0] = 1.0; u[0] = z[0] = 0.0
    for i in xrange(1,n):
        L[i] = 2*(X[i+1]-X[i-1]) - h[i-1]*u[i-1]
        u[i] = h[i]/L[i]
        z[i] = (alpha[i]-h[i-1]*z[i-1])/L[i]
    L[n] = 1.0; z[n] = c[n] = 0.0
    for j in xrange(n-1, -1, -1):
        c[j] = z[j] - u[j]*c[j+1]
        b[j] = (a[j+1]-a[j])/h[j] - (h[j]*(c[j+1]+2*c[j]))/3
        d[j] = (c[j+1]-c[j])/(3*h[j])
    splines = []
    for i in xrange(n):
        splines.append((a[i],b[i],c[i],d[i]))
    return splines,X

def plot_spline(spl,X,res=10):
    xout=[]
    yout=[]
    for i in range(len(spl)):
        x=np.linspace(X[i],X[i+1],res)
        xx=x-X[i]
        a,b,c,d=spl[i]
        y=a+b*xx+c*(xx**2)+d*(xx**3)
        xout.append(x)
        yout.append(y)
    xout=np.array(xout).flatten()
    yout=np.array(yout).flatten()
    return (xout,yout)

def newton(f,df,x,Niter=20):
    for i in range(Niter):
        x=x-f(x)/df(x)
    return x

def gaussian(sigma,norm=1., mu=0.): #returns a gaussian with height c, mean mu and width sigma as a function
    var=copy(sigma)**2
    pre=norm/np.sqrt(2*np.pi*var)
    return lambda x: pre*np.exp(-((x-mu)**2)/(2*var))

def triangle(left,right=False,norm=1.):
    if not right:
        right=left
    h=norm*2/(left+right)
    return lambda x: h*(1+x/left)*top_hat(x,-left,0)+h*(1-x/right)*top_hat(x,0,right)

