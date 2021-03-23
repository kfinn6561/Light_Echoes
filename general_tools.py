'''
Created on 14 Jan 2013

simple functions for things like progress bars and saving data to pickle files

@author: Kieran Finn
'''
import sys
import time
import pickle
import numpy as np
from copy import copy

def overprint(s):
    sys.stdout.write('\r')
    sys.stdout.flush()
    sys.stdout.write(s)

def pload(fname):
    f=open(fname,'rb')
    try:
        out=pickle.load(f)
    except:
        f.close()
        f=open(fname,'r')
        out=pickle.load(f)
    f.close()
    return out

def pdump(obj,fname):
    f=open(fname,'wb')
    pickle.dump(obj,f)
    f.close()
    
def add_comma(number): #makes a large number more readable by adding a comma every three digits
    out=''
    i=1
    number=str(number)
    while i<=len(number):
        out=number[-i]+out
        if i%3==0 and i!=len(number):
            out=','+out
        i+=1
    return out

def date_string():
    s=time.ctime()
    wday,month,day,t,year=s.split()
    t=t.split(':')
    t=t[0]+'_'+t[1]
    return day+'_'+month+'_'+year+'_'+t

def hms(time): #converts a number in seconds into a string of the form HH:MM:SS
    hours=int(time/3600)
    time-=hours*3600
    minutes=int(time/60)
    seconds=int(time-minutes*60)
    out= '%d:%02d:%02d' %(hours,minutes,seconds)
    return out
    

def progress_bar(current,total):
    screen_length=80.0
    middle=int(screen_length/2)
    timestr=' %s/%s ' %(add_comma(current+1),add_comma(total))
    timelen=len(timestr)
    start=middle-timelen/2
    end=start+timelen
    n=int(screen_length*float(current)/total)
    if n<=start:
        out='#'*n+' '*(start-n)+timestr
    elif n<=end:
        out='#'*start+timestr
    else:
        out='#'*start+timestr+'#'*(n-end)
    overprint(out)
    
def increase_sampling(x,N):
    '''takes a list samples, presumably evenly spaced,
     and returns a new list with the same range but N times the samples'''
    return np.linspace(x[0],x[-1],N*len(x))
    
def reduce_range(x,minx,maxx,include_ends=False,indices=False):
    '''creates a new list of samples restricted to the range minx-maxx.
    Both minx and maxx are included if include_ends is True'''
    i=0
    while i<len(x) and x[i]<minx:
        i+=1
    lo=copy(i)
    while i<len(x) and x[i]<maxx:
        i+=1
    hi=copy(i)
    out=copy(np.array(x[lo:hi]))
    if include_ends:
        out=np.insert(out,0,minx)
        out=np.insert(out,len(out),maxx)
    if indices:
        return [lo,hi]
    else:
        return out
    
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def extrapolate(x,y,x_new,Npoints=500):#extends a list of x,y values to x_new by linear extrapolation of the last two points
    x1,x2=x[-2:]
    y1,y2=y[-2:]
    m=(y1-y2)/(x1-x2)
    c=y1-m*x1
    
    x_new=np.linspace(x2,x_new,Npoints)[1:]
    y_new=m*x_new+c
    return(np.append(x,x_new),np.append(y,y_new))

def mat_x_vec(mat,vec):
    return np.array(mat*np.matrix(vec).T).T[0]
    
def make_steps(x,y):
    x=np.array(x)
    y=np.array(y)
    dx=x[1:]-x[:-1]
    xx=x[:-1]+dx/2.
    out_x=np.concatenate([[x[0]-dx[0]/2.],np.repeat(xx,2),[x[-1]+dx[-1]/2]])
    out_y=np.repeat(y,2)
    return (out_x,out_y)
    
    
    