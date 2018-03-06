clf
close all
clear all
clc
pkg load control

K = 1;
taua = -1;
taub = 0;
tau1 = 12;
tau2 = 3;
tau3 = 0.2;
tau4 = 0.05;
theta = 1;

%Actual function:
s = tf('s');
G = K*(taua*s+1)*(taub+1)*exp(-s)/((tau1*s+1)*(tau2*s+1)*(tau3*s+1)*(tau4*s+1));

%Taylor approximation:
thetaT = tau2 + tau3 + tau4 - taua - taub + theta
tauT = tau1
GT = K*exp(-thetaT*s)/(tauT*s+1);

%Skogestad approximation (1):
thetaS1 = tau2/2 + tau3 + tau4 - taua - taub + theta
tauS1 = tau1 + tau2/2
GS1 = K*exp(-thetaS1*s)/(tauS1*s+1);

%Skogestad approximation (2):
thetaS2 = tau3/2 + tau4 - taua - taub + theta
tauS = tau1;
tauSS = tau2 + tau3/2;
tauS2 = [tauS tauSS]
GS2 = K*exp(-thetaS2*s)/((tauS*s+1)*(tauSS*s+1));

step(G, GT, GS1, GS2)
legend('actual','Taylor','FOPDT','SOPDT','Location','SouthEast')
ylim([0 1]);
xlim([0 30]);
