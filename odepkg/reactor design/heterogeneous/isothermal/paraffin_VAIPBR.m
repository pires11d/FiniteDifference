% Variable Area Isothermal PBR

clf
close all
clear all
clc

paraffin_VAIPBR_data

[Z,Y] = ode23(@paraffin_VAIPBR_ODE,[zo zf],[Xo yo Xo yo]);

xx = Y(:,1);
yy = Y(:,2);
xxx = Y(:,3);
yyy = Y(:,4);

figure(1)
plot(Z,xx,'k-',Z,yy,'m-',Z,xxx,'k:',Z,yyy,'m:')
grid on
ylabel('X')
xlabel('z (m)')
xlim([0 zf])
ylim([0 Xf])
legend('X','y','Xpbr','ypbr','Location','SouthEast')

