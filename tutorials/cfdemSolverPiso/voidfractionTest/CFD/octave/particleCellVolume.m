%test
close all;
clear;
clc;

%====================================%
% simulation data
%====================================%
path = '../postProcessing/particleCellVolume/0/particleCellVolume.txt';
data = load(path);
Vp_sim = data(:,2);
t_sim = data(:,1);

path = '../postProcessing/volWeightedAverage/0/volWeightedAverage.txt';
data = load(path);
volAvDdtVoidFraction_sim = data(:,2);
t_sim = data(:,1);

%====================================%
dp=2*0.0015;
vp=0.3
deltaT=0.001 % coupling time
%====================================%
% analytical calculation
%====================================%
Vp=dp*dp*dp*pi/6;
np=1000;
Vptot=np*Vp

kernelSize=8
Vc=(0.1/10)^3*kernelSize
Vp=dp^3*pi/6
deltaVpdt=Vp*(deltaT*vp)/dp
avg_ddt_voidfraction=deltaVpdt/(Vc*deltaT)

%====================================%
% plot data
%====================================%
figure(1)
plot(t_sim,(Vp_sim./Vptot)*100,'r-')
hold on;
%legend("error")
title("particle volume represented in cells");
grid on;
xlabel("time in s");
ylabel("particle cell volume / particle volume in %");
print('-dpng','-r450', 'particleCellVolume')

figure(2)
plot(t_sim,volAvDdtVoidFraction_sim,'r-',[min(t_sim),max(t_sim)],[avg_ddt_voidfraction,avg_ddt_voidfraction],'b-')
grid on;
%axis([0,max(t_sim),0,100])
hold on;
legend("simulation",'analytic')
title("average ddt(voidfraction)");
grid on;
xlabel("time in s");
ylabel("average ddt(voidfraction)");
print('-dpng','-r450', 'averageDDTvoidfraction')
