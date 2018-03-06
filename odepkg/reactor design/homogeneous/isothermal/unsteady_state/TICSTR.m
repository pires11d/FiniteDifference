clf
close all
clear all
clc

TICSTR_data

to = 0;
tf = 180;

Co = [Vo Cao Cbo Cco Cdo];

[tt,ff] = ode45(@TICSTR_ODE,[to tf],Co);

%% Redefinição das Variáveis:
V = ff(:,1);
Ca = ff(:,2);
Cb = ff(:,3);
Cc = ff(:,4);
Cd = ff(:,5);

%% Equações Adicionais:
A = (pi*D^2)/4;
Q = kv*sqrt(V/A);
ma = Ca.*MMa.*Q;
mb = Cb.*MMb.*Q;
mc = Cc.*MMc.*Q;
md = Cd.*MMd.*Q;

%% Gráficos:
figure(1)
subplot(2,2,1)
    plot(tt,V)
    ylabel('V(L)')
    xlabel('t(s)')
    legend('V(t)','Location','NorthWest')
subplot(2,2,2)
    plot(tt,Ca,tt,Cb,tt,Cc,tt,Cd)
    ylabel('C(mol/L)')
    xlabel('t(s)')
    legend('Ca(t)','Cb(t)','Cc(t)','Cd(t)','Location','NorthEast')
subplot(2,2,3)
    plot(tt,Q)
    ylabel('Q(L/s)')
    xlabel('t(s)')
    legend('Q(t)','Location','NorthWest')
subplot(2,2,4)
    plot(tt,ma,tt,mb,tt,mc,tt,md)
    ylabel('m(kg/s)')
    xlabel('t(s)')
    legend('ma(t)','mb(t)','mc(t)','md(t)','Location','NorthEast')



