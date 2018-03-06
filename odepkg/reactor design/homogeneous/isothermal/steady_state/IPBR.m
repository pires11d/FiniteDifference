clf
close all
clear all
clc

IPBR_data

[V,X] = ode23(@IPBR_ODE,[Vo Vf],Xo);

figure(1)
plot(V,X,'k-')
ylabel('X')
xlabel('V (L)')
xlim([0 Vf])
ylim([0 1])
legend('X','Location','NorthEast')
