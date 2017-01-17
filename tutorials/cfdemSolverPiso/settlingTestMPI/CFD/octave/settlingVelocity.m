%
close all;
clear;
clc;

%====================================%
% simulation data
%====================================%
path = '../../DEM/post/velocity.txt';
data = load(path);
U_sim = data(:,2:4);
t_sim = data(:,1);
fprintf('final velocity of sim = %f/%f/%f m/s\n',U_sim(length(U_sim(:,1)),1),U_sim(length(U_sim(:,1)),2),U_sim(length(U_sim(:,1)),3) )


%====================================%
U_=0
X_=0
dt=0.001
tEnd=0.2
nuc = 1e-05
rhoc = 10
d_ = 0.0001
rhop = 3000
g=9.81
Uc=0
%====================================%
% analytical calculation
%====================================%
count=1;
for t=0:dt:tEnd
    count=count+1;
    t_(count)=t;
    magUr = (U_(count-1)-Uc);
    ReFunc = 1.0;
    Re = magUr*d_/nuc;
    if Re > 0.01
        ReFunc += 0.15*Re^0.687;
    end

    Dc = (24.0*nuc/d_)*ReFunc*(3.0/4.0)*(rhoc/(d_*rhop));
    U_(count) = (U_(count-1) + dt*(Dc*Uc + (1.0 - rhoc/rhop)*g))/(1.0 + dt*Dc);
    X_(count) = X_(count-1) + dt*U_(count);
    Re_(count) = Re; 

end
fprintf('final velociy = %f m/s\n',U_(length(U_)))
fprintf('final position = %f m\n',X_(length(X_)))

%====================================%
% plot data
%====================================%
figure(1)
plot(t_,U_, 'k--')
hold on
plot(t_sim,-U_sim(:,2),'rd-')
legend("analytical - Stokes","simulation - DiFelice?") 

%print('cfdemSolverPiso_settlingTestMPI.eps','-deps2')
print -color "cfdemSolverPiso_settlingTestMPI.png"

%figure(2)
%plot(t_,X_)

%figure(3)
%plot(t_,Re_)
