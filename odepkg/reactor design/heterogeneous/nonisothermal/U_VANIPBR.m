% Variable Area Nonisothermal PBR

clf
close all
clear all
clc

U_VANIPBR_data

[Z,Y] = ode23s(@U_VANIPBR_ODE,[zo zf],[Xo Xo yo yo To To Wo Wo]);

xx = Y(:,1);
xxx = Y(:,2);
yy = Y(:,3);
yyy = Y(:,4);
TT = Y(:,5);
TTT = Y(:,6);
W = Y(:,7);
W2 = Y(:,8);
Pa = Pao*((1-xx)/(1+e*xx))*yy;
Pa2 = Pao*((1-xxx)/(1+e*xxx))*yyy;

figure(1)
plot(Z,xx,'k-',Z,yy,'m-',Z,xxx,'k:',Z,yyy,'m:')
grid on
ylabel('X')
xlabel('z (m)')
xlim([0 zf])
ylim([0 Xf])
legend('X','y','X.pbr','y.pbr','Location','SouthEast')

figure(2)
plot(Z,TT,'r-',Z,TTT,'r:')
grid on
ylabel('T (K)')
xlabel('z (m)')
xlim([0 zf])
ylim([To-10 To+10])
legend('T','T.pbr','Location','NorthEast')

figure(3)
plot(Z,W,'b-',Z,W2,'b:')
grid on
ylabel('W (kg)')
xlabel('z (m)')
xlim([0 zf])
ylim([0 max(W)])
legend('W','W.pbr','Location','NorthWest')

figure(4)
plot(Z,Pa,'k-',Z,Pa2,'k:')
grid on
ylabel('Pa (bar)')
xlabel('z (m)')
xlim([0 zf])
ylim([0 Po])
legend('Pa','Pa.pbr','Location','SouthWest')



