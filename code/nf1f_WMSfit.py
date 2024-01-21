#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 14:51:03 2019

@author: saisu
""",
import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter
import scipy.signal as signal

def nf1f_WMSfit(v,*args): 
    x=args
    #global num
    #num=num+1
    wV=0.5346*x[1]+(0.2166*x[1]**2+x[0]**2)**0.5
    drt_ved=x[0]/2/(m.sqrt(m.log(2)))
    beta=drt_ved/(x[1]/2+drt_ved)  
    faiv=(beta/drt_ved/m.sqrt(m.pi))+2*((1-beta)/(m.pi*x[1]))
    t=np.linspace(0,1,len(v))
    fm=350            #modulation frequency
    fai_0=0.5*m.pi     #phase shift of FM
    a=x[4]              #modulation depth
    temp=2*m.pi*fm*t+fai_0
    temp=[m.cos(x) for x in temp]
    v_mod=v+np.dot(a,temp) 
    tempt=(np.multiply(-4*m.log(2),((v_mod-x[5])/wV)**2))
    tempt2=-0.4*abs(abs((v_mod-x[5])/wV)**2.25)
    tempt=[m.exp(x) for x in tempt]
    tempt2=[m.exp(x) for x in tempt2]
    vv=(np.dot((1-x[1]/wV),tempt)+\
    x[1]/wV/(4*((v_mod-x[5])/wV)**2+1)+\
    0.016*(1-x[1]/wV)*x[1]/wV *(tempt2-10/(abs(abs(v_mod-x[5]/wV)**2.25)+10)))*faiv
    absorb=x[2]*vv

    #I0 after modulation
    I_max=0.9
    I_min=0
    k=(I_max-I_min)/1.2
    I0_aver=I_min+np.dot(k,(t+0.2)); #ramp base
    i0=0.1              #amplitude of linear IM
    fai_1=x[3]          #phase shift of FM/IM
    i2=0                #amplitude of nonlinear IM
    fai_2=1.31*m.pi     #phase shift of nonlinear IM
    temp=2*m.pi*fm*t+fai_1
    temp=[m.cos(x) for x in temp]
    temp2=4*m.pi*fm*t+fai_2
    temp2=[m.cos(x) for x in temp2]
    I0=I0_aver*(1+i0/I0_aver*temp+i2/I0_aver*temp2);
    
    absorb=np.multiply(-1,absorb)
    tempt=[m.exp(x) for x in absorb]
    It=np.multiply(I0,tempt) #transimitted laser intensity

    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t,It)
    '''

    #the following program are the lock-in
    sig=It
    n=np.linspace(1,len(t),len(t))
    N=len(t)/(max(t)-min(t))/fm
    
    # 2f
    harmonic=2
    temp=harmonic*2*m.pi*n/N
    Vrs=[m.sin(x) for x in temp]
    Vrc=[m.cos(x) for x in temp]
    x2=np.multiply(sig,Vrs)
    y2=np.multiply(sig,Vrc)
    b, a = signal.butter(4,1/7.6/N, 'low')
    X2 = lfilter(b,a, x2)
    Y2 = lfilter(b,a, y2)
    tempt=X2**2+Y2**2
    tempt=[m.sqrt(x) for x in tempt]
    S2=np.dot(2,tempt)
    
    # 1f
    harmonic=1
    temp=harmonic*2*m.pi*n/N
    Vrs=[m.sin(x) for x in temp]
    Vrc=[m.cos(x) for x in temp]
    x1=np.multiply(sig,Vrs)
    y1=np.multiply(sig,Vrc)
    X1 = lfilter(b,a, x1)
    Y1 = lfilter(b,a, y1)
    tempt=X1**2+Y1**2
    tempt=[m.sqrt(x) for x in tempt]
    S1=np.dot(2,tempt)
    norm_2f_calculated=[S2[k]/S1[k] for k in range(len(S2))] 

    return norm_2f_calculated[int(len(t)/5): int(len(t))-int(len(t)/5)]












