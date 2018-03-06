% Nonisothermal PBR

clear
clc
clf

nitrous_data

[Z,Y] = ode23s(@nitrous_ODE,[zo zf],[Xo Xo yo yo To To Wo Wo]);

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
Pb = Pao*((thetab-(b/a)*xx)/(1+e*xx))*yy;
Pb2 = Pao*((thetab-(b/a)*xxx)/(1+e*xxx))*yyy;

figure(1)
subplot(2,2,1)
plot(Z,xx,'k-',Z,yy,'m-',Z,xxx,'k:',Z,yyy,'m:')
grid on
ylabel('X, y')
xlabel('z (m)')
xlim([0 zf])
ylim([0 Xf])
legend('X','y','X.pbr','y.pbr','Location','SouthEast')

subplot(2,2,2)
plot(Z,TT,'r-',Z,TTT,'r:')
grid on
ylabel('T (K)')
xlabel('z (m)')
xlim([0 zf])
ylim([To-10 To+50])
legend('T','T.pbr','Location','NorthEast')

subplot(2,2,3)
plot(Z,W,'b-',Z,W2,'b:')
grid on
ylabel('W (kg)')
xlabel('z (m)')
xlim([0 zf])
ylim([0 max(W2)])
legend('W','W.pbr','Location','NorthWest')

subplot(2,2,4)
plot(Z,Pa,'k-',Z,Pa2,'k:',Z,Pb,'g-',Z,Pb2,'g:')
grid on
ylabel('Pn2o,Pno (bar)')
xlabel('z (m)')
xlim([0 zf])
ylim([0 Pbo])
legend('Pn2o','Pn2o.pbr','Pno','Pno.pbr','Location','NorthEast')



