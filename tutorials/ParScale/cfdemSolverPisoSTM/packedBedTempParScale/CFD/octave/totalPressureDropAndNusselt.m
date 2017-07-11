close all;
clear;
clc;


%*********************************************************************%
% pressure drop
%*********************************************************************%

%===================
% Ergun Equation
%===================
fprintf('\ncalc Ergun eqn:\n')
dp = 0.022;			% particle diameter
phip = 1;			% sphericity
epsilon = 0.4436;		        % voidfraction
U = 1;    			% superficial velocity
L = 1;			    % length of bed
rhoG = 1.188;		% density in kg/m
nuG = 1.5e-3;	    % kinemat Visk in m2/s
muG = nuG*rhoG		% dynam visc in Pa s

dpErgun= L * (
                150*((1-epsilon)^2/epsilon^3)*((muG*U)/(phip*dp)^2) 
              +1.75*((1-epsilon)/epsilon^3)*((rhoG*U^2)/(phip*dp))
        )

fprintf('final pressure drop = %f Pa\n',dpErgun)
%====================================%
% simulation data
%====================================%
%path = '../probes/0/p'; % ext32
path = '../postProcessing/probes/0/p';

data = load(path);
[x,y]=size(data);
dp_sim = data(:,2)-data(:,y);
t_sim = data(:,1);
fprintf('final pressureDrop of sim = %f Pa\n',dp_sim(length(dp_sim)) )

%====================================%
% plot data
%====================================%
figure(1)
plot(t_sim,dpErgun*ones(1,length(t_sim)),t_sim,dp_sim)
axis([0,0.7,0,dpErgun(length(dpErgun))])
title("Ergun pressure drop")
legend("analytical - Ergun","simulation") 

%print('cfdemSolverPiso_settlingTest.eps','-deps2')
print -color "cfdemSolverPisoSTM_pressureDrop.png"

%*********************************************************************%
% heat transfer
%*********************************************************************%

%====================================%
% simulation data
%====================================%

cp = 1007;
A = 0.01;
U = U;
alpha = epsilon;
nuF = nuG;
rhoF = rhoG;
Tp = 600;
Np = 1005;
lambda = 0.0256;

%path = '../probes/0/T'; % ext32
path = '../postProcessing/probes/0/T';

data = load(path);
[x,y]=size(data);
Tin_sim =  data(:,2);                                          % mean temp inlet temp [K]
Tout_sim = data(:,3);                                          % mean temp outlet temp [K]
t_sim = data(:,1);
fprintf('final inlet temperature of sim = %f K\n',Tin_sim(length(Tin_sim)) )
fprintf('final outlet temperature of sim = %f K\n',Tout_sim(length(Tout_sim)) )


uF = U/alpha;                                               % interstitial velocity [m/s]
ReP_sim = uF.*ones(length(t_sim),1)*dp/nuF;                                       % ReynoldsNr based in dp
Pr = nuF*rhoF*cp/lambda;
qin_sim = U * A * rhoF * cp .* Tin_sim;
qout_sim = U * A * rhoF * cp .* Tout_sim;
q_sim = (qout_sim-qin_sim);                                 % particle fluid heat flux [W] (out-in)
Tmean_sim = 0.5*(Tin_sim+Tout_sim);                               % mean temp of fluid
deltaT = Tp - Tmean_sim;                                          % mean temp diff between partcles and fluid
h=q_sim./(Np*dp^2*pi*deltaT);                                      % average particle-fluid heat transfer coeff [W/(m2*K)]
Nu_sim = h.*dp/lambda;                                             % mean particle Nusselt nr
t_sim = data(:,1);
fprintf('q_sim = %f \n',q_sim(length(Nu_sim)))
fprintf('Nu_sim = %f \n',Nu_sim(length(Nu_sim)))

%===================
% Nusselt Equation
%===================
fprintf('\ncalc Nusselt eqn:\n')

%% following Li and Mason
n=3.5
if (ReP_sim(length(uF)) <200)
     Nu_LiMason=2.*ones(length(uF),1)+0.6*alpha^n.*ReP_sim.^0.5.*Pr^0.33;
elseif (ReP_sim(length(uF)) <1500)
     Nu_LiMason=2.*ones(length(uF),1)+0.5*alpha^n.*ReP_sim.^0.5*Pr^0.33 + 0.02*alpha^n.*ReP_sim.^0.8.*Pr^0.33;
else
     Nu_LiMason=2.*ones(length(uF),1) +0.000045*alpha^n.*ReP_sim.^1.8;
end

fprintf('Nu_LiMason = %f \n',Nu_LiMason(length(Nu_LiMason)))

%====================================%
% plot data
%====================================%
figure(2)
plot(t_sim,Nu_LiMason,t_sim,Nu_sim)
title("Nusselt nr")
legend("analytical - ","simulation") 

%print('cfdemSolverPisoScalar_NusseltNr.eps','-deps2')
print -color "cfdemSolverPisoSTM_Nusselt.png"

figure(3)
plot(t_sim,Tin_sim,t_sim,Tout_sim)
title("inlet/outlet temperature")
legend("inlet","outlet") 

%print('cfdemSolverPisoScalar_NusseltNr.eps','-deps2')
print -color "cfdemSolverPisoSTM_temperatures.png"

