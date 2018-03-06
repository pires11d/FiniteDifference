clf
close all
clear all
clc

pkg load odepkg

reaction2_data

C0 = [dCao Cao dCbo Cbo dCro Cro];

[x,C] = ode23(@reaction2_ode,[0 xf],C0);

%% Re-Definitions
dCa = C(:,1);
Ca = C(:,2);
dCb = C(:,3);
Cb = C(:,4);
dCr = C(:,5);
Cr = C(:,6);

%% Plots:
figure(1)
    plot(x,Ca,x,Cb,x,Cr)
    grid on
    ylabel('C (mol/L)')
    xlabel('x (cm)')
    legend('Ca(x)','Cb (x)','Cr (x)')
    xlim([0 xf])
    ylim([0 1])
