#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:02:54 2019
@author: saisu
"""
from scipy.io import loadmat
import math as m
import numpy as np
import matplotlib.pyplot as plt
import nf1f_WMSfit
import scipy.signal as signal
from scipy.signal import lfilter
import scipy.optimize as optimization
import time
import warnings

warnings.filterwarnings("ignore")

def downsampling(data,space):
    new_list=[]
    k=int(len(data)/space)
    for i in range(k):
        new_list.append(data[i*int(space)])
    return new_list
    

################################# Main #############################
num=0
H2O_database= loadmat("/Users/saisu/Desktop/Final Year Projrct/H2O_300K.mat")
M=18 #molecule H2O
P=1.00 #pressure 1 atm
T=300 #temperature K (must fixed here)
x_s=0.02 #mole fraction of H2O
x_air=1-x_s # mole fraction of air
L=50 #pathlength cm
v0=7407.8  #center frequency cm-1
res_v=0.0002  #spectra resolution in calculation;
v0_L=1.00
v0_R=1-res_v    #selected spectra range in calculation; v0_L: range left to v0;  v0_R: range right to v0
value_list=list(H2O_database.values())
H2O_300K=value_list[3] 
col1_v=H2O_300K[:,0].tolist(); col2_St=H2O_300K[:,1].tolist();gamma_air=H2O_300K[:,2].tolist();gamma_self=H2O_300K[:,3].tolist();col5_E=H2O_300K[:,4].tolist();

index_v=list()  # initiate empty list
index_St=list()
index_gamma_air=list()
index_gamma_self=list()
index_E=list()
wD=list()
wC=list()
wV=list()
drt_ved=list()
beta=list()
faiv0=list()
#for i in H2O_300K[:,0]: #iterate all v, i.e col1:v
for i in range(len(col1_v)):
   if col1_v[i]>=v0-v0_L and col1_v[i]<=v0+v0_R:
     index_v.append(col1_v[i])  # find the transitions in this range
     index_St.append(col2_St[i]) 
     index_gamma_air.append(gamma_air[i]) 
     index_gamma_self.append(gamma_self[i])
     index_E.append(col5_E[i]) 
     
##############     calculate the parameters of the selected transitions     ###############
#for i in range(len(col1_v)):
for i in range(len(index_v)):
     wD.append(np.multiply(7.1623e-7,index_v[i]*(T/M)**0.5))        #Dopple width
     wC.append(P*2*(x_s*index_gamma_self[i]+x_air*index_gamma_air[i]))      #Collision width
     wV.append(0.5346*wC[i]+(0.2166*wC[i]**2+wD[i]**2)**0.5)                 #Voigt width
     drt_ved.append(wD[i]/2/m.sqrt(m.log(2)))
     beta.append(drt_ved[i]/(wC[i]/2+drt_ved[i]))
     faiv0.append((beta[i]/drt_ved[i]/m.sqrt(m.pi))+2*((1-beta[i])/(m.pi*wC[i])))
 
# v before modulation
v=np.linspace(v0-v0_L,v0+v0_R,10000)
# v after modulation
t=np.linspace(0,1,len(v))
fm=350          #modulation  frequency
N=len(v)/fm     #number of samples per cycle
fai_0=0.5*m.pi  #phase shift of FM
a=0.1           #modulation depth
temp=2*m.pi*fm*t+fai_0
temp=[m.cos(x) for x in temp]
v_mod=v+np.dot(a,temp) 

'''
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t,v_mod)
plt.xlabel('Time/s')
plt.ylabel('Frequency/cm-1')
plt.title('Frequency Response to the Laser Injection Current')
plt.show()
'''

###################    DAS voigt lineshape     ##########################
vv_DAS_per_v=np.mat(np.zeros((len(index_v),len(v))))

for j in range(len(index_v)): #20  
    temp1=-4*m.log(2)*((v-index_v[j])/wV[j])**2
    temp1=[m.exp(x) for x in temp1]
    temp2=-0.4*abs(abs((v-index_v[j])/wV[j])**2.25)
    temp2=[m.exp(x) for x in temp2]
    vv_DAS_per_v[j,:]=(P*x_s*L*index_St[j]*7.34e21/T*(np.multiply((1-wC[j]/wV[j]),temp1)+wC[j]/wV[j]/(4*((v-index_v[j])/wV[j])**2+1)+0.016*(1-wC[j]/wV[j])*wC[j]/wV[j]*(temp2-10/(abs(abs((v-index_v[j])/wV[j])**2.25)+10)))*faiv0[j])
    damn=np.trapz(vv_DAS_per_v[j,:],v)
    vv_DAS_per_v[j,:]=np.multiply(vv_DAS_per_v[j,:],P*x_s*L*index_St[j]*7.34e21/T/np.trapz(vv_DAS_per_v[j,:],v))
absorb_DAS=np.array( vv_DAS_per_v.sum(axis=0))

        
#Figure Plot    
'''
plt.figure()
plt.plot(v,absorb_DAS,label="Simulation (Voigt)")
plt.legend(loc='upper left')
plt.text(7406.85, 0.014, "Molecule: H2O", size = 12, alpha = 0.9)
plt.text(7406.85, 0.012, "T: 300K", size = 12, alpha = 0.9)
plt.text(7406.85, 0.010, "L: 50cm", size = 12, alpha = 0.9)
plt.text(7406.85, 0.008, "P: 1atm", size = 12, alpha = 0.9)
plt.text(7406.85, 0.006, "v0: 7407.8cm-1", size = 12, alpha = 0.9)
plt.text(7406.85, 0.004, "X: 2%", size = 12, alpha = 0.9)
plt.xlabel('Frequency/cm-1')
plt.ylabel('Absorbance')
plt.title('Absorbance vs. Frequency')
plt.show()
'''

#################       Modulated DAS voigt lineshape      #############################

vv_DAS_modulated=np.mat(np.zeros((len(index_v),len(v_mod))))

for j in range(len(index_v)): #20  
    temp1=-4*m.log(2)*((v_mod-index_v[j])/wV[j])**2
    temp1=[m.exp(x) for x in temp1]
    temp2=-0.4*abs(abs((v_mod-index_v[j])/wV[j])**2.25)
    temp2=[m.exp(x) for x in temp2]
    vv_DAS_modulated[j,:]=(P*x_s*L*index_St[j]*7.34e21/T*(np.multiply((1-wC[j]/wV[j]),temp1)+wC[j]/wV[j]/(4*((v_mod-index_v[j])/wV[j])**2+1)+0.016*(1-wC[j]/wV[j])*wC[j]/wV[j]*(temp2-10/(abs(abs((v_mod-index_v[j])/wV[j])**2.25)+10)))*faiv0[j])
    vv_DAS_modulated[j,:]=np.multiply(vv_DAS_modulated[j,:],P*x_s*L*index_St[j]*7.34e21/T/np.trapz(vv_DAS_modulated[j,:],v_mod))
absorb_modulated=np.array(vv_DAS_modulated.sum(axis=0)).T


#Figure Plot    
'''
plt.figure()
plt.plot(t,absorb_modulated,label="Simulation")
plt.legend(loc='upper left')
plt.xlabel('Time/s')
plt.ylabel('Absorbance')
plt.title('Modulated Absorbance vs. Time')
plt.show()
'''
        
#######################  I0 after modulation   #########################
I_max=0.9
I_min=0
k=(I_max-I_min)/1.2
I0_aver=I_min+np.dot(k,(t+0.2))  #ramp base
i0=0.1              #amplitude of linear IM
fai_1=-0.2*m.pi          #phase shift of FM/IM
i2=0                #amplitude of nonlinear IM
fai_2=1.31*m.pi     #phase shift of nonlinear IM
temp=2*m.pi*fm*t+fai_1
temp=[m.cos(x) for x in temp]
temp2=4*m.pi*fm*t+fai_2
temp2=[m.cos(x) for x in temp2]
I0=I0_aver*(1+i0/I0_aver*temp+i2/I0_aver*temp2)
absorb_modulated=np.multiply(-1,absorb_modulated)
tempt=[m.exp(x) for x in absorb_modulated]
It=np.multiply(I0,tempt) #transimitted laser intensity
 #Figure Plot    

plt.figure()
plt.plot(t,It,label="Simulation")  
plt.title('Absorption-free Laser Intensity vs. time')     


###############   Lock-in detection    #################################
sig=It
n=np.linspace(1,len(t),len(t))

# 2f
harmonic=2
temp=harmonic*2*m.pi*n/N
Vrs=[m.sin(x) for x in temp]
Vrc=[m.cos(x) for x in temp]
x2=np.multiply(sig,Vrs)
y2=np.multiply(sig,Vrc)
#b, a = signal.butter(4,1/7.6/N, 'low')
order, wn = signal.buttord(1/8/N, 1/2/N,gpass=1/3.5,gstop=60, analog='true')
b, a = signal.butter(order,wn, 'low')
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
norm_2f=[S2[k]/S1[k] for k in range(len(S2))]

'''
plt.figure()
plt.plot(v[int(len(t)/5): int(len(t))],norm_2f[int(len(t)/5): int(len(t))],label="Simulation")  
plt.legend(loc='upper left')
plt.xlabel('Frequency/ cm-1')
plt.ylabel('S2f/1f')
plt.title('1f-noralized 2f signal  vs. frequency')     
'''

###############    fitting the CF-WMS 2f/1f    ###############################

x0=[0.05, 0.1, 0.001, 0, 0.05,v0]          #Starting guess
lb=[0.01, 0.08, 0.0001, -1*m.pi, 0.05, v0-v0_L]    #lower bound
ub=[0.1, 0.15, 0.1, m.pi, 0.2, v0+v0_R]     #upper bound
    #   x is the vector of the fitting parameters, x(0)=wD; x(1)=wC;
    #   x(2)=P*X*S(T)*L; x(3)=fai_1; x(4)=a; x(5)=v0
bdy=[lb,ub] # boundary

input_data=v
output_data=norm_2f[int(len(t)/5): int(len(t))-int(len(t)/5)]

new_in=downsampling(input_data, 4)
new_out=downsampling(output_data,4)

start_3 = time.clock()
popt,pconv=optimization.curve_fit(nf1f_WMSfit.nf1f_WMSfit,new_in,new_out,x0,bounds=bdy,ftol=1e-1)
x=popt
elapsed_3 = (time.clock() - start_3)
print("Time used for Fitting:",elapsed_3)

###############    fitting the norm_2f_1f with fitted parameters    ###############################
norm_2f_fitted=nf1f_WMSfit.nf1f_WMSfit(v,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])

###############    fitting the absorbance with fitted parameters    ###############################

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
vv_mod=(np.dot((1-x[1]/wV),tempt)+\
x[1]/wV/(4*((v_mod-x[5])/wV)**2+1)+\
0.016*(1-x[1]/wV)*x[1]/wV *(tempt2-10/(abs(abs(v_mod-x[5]/wV)**2.25)+10)))*faiv
absorb_modulated=x[2]*vv_mod    #modulated absorbance

tempt=(np.multiply(-4*m.log(2),((v-x[5])/wV)**2))
tempt2=-0.4*abs(abs((v-x[5])/wV)**2.25)
tempt=[m.exp(x) for x in tempt]
tempt2=[m.exp(x) for x in tempt2]
vv=(np.dot((1-x[1]/wV),tempt)+\
x[1]/wV/(4*((v-x[5])/wV)**2+1)+\
0.016*(1-x[1]/wV)*x[1]/wV *(tempt2-10/(abs(abs(v-x[5]/wV)**2.25)+10)))*faiv
absorb_fit=x[2]*vv              #fitted absorbance



#Figure Plot 
plt.figure()
plt.plot(norm_2f[int(len(t)/5): int(len(t))-int(len(t)/5)],color='r',label="Measured")
plt.plot(norm_2f_fitted,color='b',label="Fitted")
plt.legend(loc='upper left')
plt.xlabel('Frequency/ cm-1')
plt.ylabel('S2f/1f')
plt.title('Fitted S2f/1f vs. Measured S2f/1f')  

    
    

    
    
    
    
    
    
    
    
    







