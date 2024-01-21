% WMS for selected H2O spectra 
close all;
clear all; clc;
load H2O_300K.mat; %col1:v; col2:S(300k); col3:gamma_air; col4:gamma_self; col5:E"
M=18; %molecule H2O
P=1; %pressure 1 atm
T=300; %temperature K (must fixed here)
x_s=0.02; % mole fraction of H2O
x_air=1-x_s; % mole fraction of air
L=50; % pathlength cm
global num;
num=0;
global damn;
damn=0;
%% plot modulated v(t) 
v0=7407.8; %center frequency cm-1
res_v=0.0002; %spectra resolution in calculation;
v0_L=1; v0_R=1-res_v; %selected spectra range in calculation; v0_L: range left to v0;  v0_R: range right to v0
index_v=find((H2O_300K(:,1)>=v0-v0_L)&(H2O_300K(:,1)<=v0+v0_R)); % find the transitions in this range

%calculate the parameters of the selected transitions
wD=(7.1623e-7)*H2O_300K(index_v,1)*(T/M).^0.5; % Dopple width
wC=P*2*(x_s*H2O_300K(index_v,4)+x_air*H2O_300K(index_v,3)); %Collision width
wV=0.5346*wC+(0.2166*wC.^2+wD.^2).^0.5; %Voigt width
drt_ved=wD/2/sqrt(log(2));
beta=drt_ved./(wC/2+drt_ved);
faiv0=(beta./drt_ved/sqrt(pi))+2*((1-beta)./(pi*wC));   % Peak height of the Voigt Lineshape

%v before modulation
v=v0-v0_L:res_v:v0+v0_R;

%v after modulation
t=linspace(0,1,length(v));
fm=200 %modulation  frequency
N=length(v)/fm %number of samples per cycle
fai_0=0.5*pi;  %phase shift of FM
a=0.1; %modulation depth
v_mod=v+a*cos(2*pi*fm*t+fai_0);
% plot(t,v_mod)

% DAS voigt lineshape
for i=1:length(index_v)
    vv_DAS(:,i)=P*x_s*L*H2O_300K(index_v(i),2)*7.34e21/T*((1-wC(i)./wV(i)).*exp(-4*log(2)*((v-H2O_300K(index_v(i),1))./wV(i)).^2)+wC(i)./wV(i)./(4*((v-H2O_300K(index_v(i),1))./wV(i)).^2+1)+0.016*(1-wC(i)./wV(i))*wC(i)./wV(i).*(exp(-0.4*abs(((v-H2O_300K(index_v(i),1))./wV(i)).^2.25))-10./(abs(((v-H2O_300K(index_v(i),1))./wV(i)).^2.25)+10)))*faiv0(i);
    vv_DAS(:,i)=vv_DAS(:,i)*P*x_s*L*H2O_300K(index_v(i),2)*7.34e21/T/trapz(v,vv_DAS(:,i));
end
absorb_DAS=sum(vv_DAS,2);
 %plot(v,absorb_DAS); hold on;
% xlabel('Wavenumber [cm-1]');ylabel('Absorbance')
A_DAS=trapz(v,absorb_DAS);

% modulated voigt lineshape
for i=1:length(index_v)
    vv(:,i)=P*x_s*L*H2O_300K(index_v(i),2)*7.34e21/T*((1-wC(i)./wV(i)).*exp(-4*log(2)*((v_mod-H2O_300K(index_v(i),1))./wV(i)).^2)+wC(i)./wV(i)./(4*((v_mod-H2O_300K(index_v(i),1))./wV(i)).^2+1)+0.016*(1-wC(i)./wV(i))*wC(i)./wV(i).*(exp(-0.4*abs(abs((v_mod-H2O_300K(index_v(i),1))./wV(i)).^2.25))-10./(abs(abs((v_mod-H2O_300K(index_v(i),1))./wV(i)).^2.25)+10)))*faiv0(i);
    vv(:,i)=vv(:,i)*P*x_s*L*H2O_300K(index_v(i),2)*7.34e21/T/trapz(v,vv(:,i));
end
absorb=sum(vv,2);
% plot(t,absorb); 

%I0 after modulation
I_max=0.9;
I_min=0;
k=(I_max-I_min)/1.2;
I0_aver=I_min+k.*(t+0.2); %ramp base
% I0_aver=0.4*cos(2*pi*1/2*t)+0.5; %sine base
i0=0.1; %amplitude of linear IM
fai_1=-0.2*pi; %phase shift of FM/IM
i2=0;%amplitude of nonlinear IM
fai_2=1.31*pi;%phase shift of nonlinear IM
I0=I0_aver.*(1+i0./I0_aver.*cos(2*pi*fm*t+fai_1)+i2./I0_aver.*cos(4*pi*fm*t+fai_2));
% figure(1); plot(t,I0); hold on;

It=I0.*exp(-absorb'); %transimitted laser intensity
% It=awgn(It,60,'measured'); %add Guassian white noise to the transmission
figure(2); plot(t,It);
 xlabel('Time');ylabel('Transmission')

%% the following program are the lock-in
sig=It;
xx=1:length(t);
harmonic=2;
Vrs=sin(harmonic*2*pi.*xx/N);
Vrc=cos(harmonic*2*pi.*xx/N);
x2=sig.*Vrs;
y2=sig.*Vrc;
d=fdesign.lowpass('Fp,Fst,Ap,Ast',1/4/N,1/N,1,60);
Hd = design(d,'butter');
X2= filter(Hd,x2);
Y2=filter(Hd,y2);
S2=2*sqrt(X2.^2+Y2.^2);
%figure();plot(t,S2);hold on;

harmonic=1;
Vrs=sin(harmonic*2*pi.*xx/N);
Vrc=cos(harmonic*2*pi.*xx/N);
x1=sig.*Vrs;
y1=sig.*Vrc;
X1= filter(Hd,x1);
Y1=filter(Hd,y1);
S1=2*sqrt(X1.^2+Y1.^2);
% plot(t,S1);

% freq=0:(2*pi)/length(x1):pi;
% x1dft=fft(x2);
% X1dft=fft(X2);
% figure
% plot(freq/pi,abs(x1dft(1:length(x1)/2+1)))
% hold on
% plot(freq/pi,abs(X1dft(1:length(X1)/2+1)))

norm_2f=S2./S1;
%figure()
%plot(t,norm_2f);
%5xlabel('Time');ylabel('S 2f/1f')

%% fitting the CF-WMS 2f/1f
x0=[0.05; 0.1;0.001;0;0.05; v0]; % Starting guess
lb=[0.01; 0.08;1e-4;-pi;0.05; v0-v0_L]; %lower bound
ub=[0.1; 0.15;0.1;pi;0.2; v0+v0_R]; %upper bound
%x is the vector of the fitting parameters, x(1)=wD; x(2)=wC;
%x(3)=P*X*S(T)*L; x(4)=fai_1; x(5)=a; x(6)=v0
options = optimoptions(@lsqcurvefit,'MaxIter',1500);
tic;
[x,resnorm,residual] = lsqcurvefit(@nf1f_WMSfit,x0,v,norm_2f((length(t)/5):end),lb,ub,options); 
toc;
wV_fit=0.5346*x(2)+(0.2166*x(2)^2+x(1)^2)^0.5;
drt_ved_fit=x(1)/2/sqrt(log(2));
beta_fit=drt_ved_fit./(x(2)/2+drt_ved_fit);
faiv_fit=(beta_fit./drt_ved_fit/sqrt(pi))+2*((1-beta_fit)./(pi*x(2)));
v_mod_fit=v+x(5)*cos(2*pi*fm*t+fai_0);
vv_fit=((1-x(2)/wV_fit)*exp(-4*log(2)*((v_mod_fit-x(6))/wV_fit).^2)+x(2)./wV_fit./(4*((v_mod_fit-x(6))./wV_fit).^2+1)+0.016*(1-x(2)/wV_fit)*x(2)/wV_fit.*(exp(-0.4*abs(((v_mod_fit-x(6))/wV_fit).^2.25))-10./(abs(((v_mod_fit-x(6))/wV_fit).^2.25)+10)))*faiv_fit;
absorb_fit=x(3)*vv_fit;
a_fit=x(3)*((1-x(2)/wV_fit)*exp(-4*log(2)*((v-x(6))/wV_fit).^2)+x(2)./wV_fit./(4*((v-x(6))./wV_fit).^2+1)+0.016*(1-x(2)/wV_fit)*x(2)/wV_fit.*(exp(-0.4*abs(((v-x(6))/wV_fit).^2.25))-10./(abs(((v-x(6))/wV_fit).^2.25)+10)))*faiv_fit;%calculate fitted lineshape
A_fit=trapz(v,a_fit);%calculate fitted integrated absorbance

figure();
plot(v,absorb_DAS,'color','r'); hold on;
xlabel('Wavenumber [cm-1]');ylabel('Absorbance')
plot(v,a_fit,'color','b');
xlabel('Wavenumber [cm-1]');ylabel('Absorbance')

I0_fit=I0_aver.*(1+i0./I0_aver.*cos(2*pi*fm*t+x(4))+i2./I0_aver.*cos(4*pi*fm*t+fai_2));
It_fit=I0_fit.*exp(-absorb_fit);
sig_fit=It_fit;

harmonic=2;
Vrs=sin(harmonic*2*pi.*xx/N);
Vrc=cos(harmonic*2*pi.*xx/N);
x2_fit=sig_fit.*Vrs;
y2_fit=sig_fit.*Vrc;
d_fit=fdesign.lowpass('Fp,Fst,Ap,Ast',1/4/N,1/N,1,60);
Hd_fit = design(d_fit,'butter');
X2_fit= filter(Hd_fit,x2_fit);
Y2_fit=filter(Hd_fit,y2_fit);
S2_fit=2*sqrt(X2_fit.^2+Y2_fit.^2);

harmonic=1;
Vrs_fit=sin(harmonic*2*pi.*xx/N);
Vrc_fit=cos(harmonic*2*pi.*xx/N);
x1_fit=sig.*Vrs_fit;
y1_fit=sig.*Vrc_fit;
X1_fit= filter(Hd_fit,x1_fit);
Y1_fit=filter(Hd_fit,y1_fit);
S1_fit=2*sqrt(X1_fit.^2+Y1_fit.^2);
 %plot(S1,'r');
 norm_2f_fit=S2_fit./S1_fit;

figure()
plot(t((length(t)/5):end),norm_2f((length(t)/5):end),'color','red');hold on;
plot(t((length(t)/5):end),norm_2f_fit((length(t)/5):end),'color','blue');%fitted WMS


noob=vv_fit./a_fit;
