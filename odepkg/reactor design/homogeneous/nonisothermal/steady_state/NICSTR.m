% Nonisothermal CSTR

clf
close all
clear all
clc

NICSTR_data

X = [0:0.05:Xf];

% Energy Balance
T = (U*A*Ta+Fao*((-dHr+dCp*Tr)*X+sumCpi*To))./(U*A+Fao*(sumCpi+dCp*X));

% Reaction Rate
k = ko*exp((Ea/R)*((1/To)-(1./T)));
Ca = Cao*((1-X)./(1+e*X)).*(To./T);
ra = k.*Ca;

% Mole Balance
V = Fao*X./ra;

figure(1)
plot(X,V,'k-')
grid on
ylabel('V (L)')
xlabel('X')
xlim([0 Xf])
ylim([0 1000])

figure(2)
plot(X,T,'k-')
grid on
ylabel('T (K)')
xlabel('X')
xlim([0 Xf])
ylim([360 380])

Tf=T(19)
Vf=V(19)

