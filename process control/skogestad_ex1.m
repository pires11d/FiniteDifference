clf
close all
clear all
clc
pkg load control

K = 1;
taua = -0.1;
taub = 0;
tau1 = 5;
tau2 = 3;
tau3 = 0.5;

%Actual function:
num = K*[taua*taub taua+taub 1];
den = [tau1*tau2*tau3 tau1*tau2+tau1*tau3+tau2*tau3 tau1+tau2+tau3 1];
G = tf(num,den)
s = tf('s');

%Taylor approximation:
thetaT = tau2 + tau3 - taua
tauT = tau1
GT = K*exp(-thetaT*s)/(tauT*s+1);

%Skogestad approximation:
thetaS = tau2/2 + tau3 - taua
tauS = tau1 + tau2/2
GS = K*exp(-thetaS*s)/(tauS*s+1);

step(G, GT, GS)
legend('actual','Taylor','FOPDT','Location','SouthEast')
ylim([0 1]);
