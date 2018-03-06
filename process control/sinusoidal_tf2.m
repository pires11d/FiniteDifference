clf
close all
clear all
clc
pkg load control

K = 1;
taua = 0;
tau1 = 1;
tau2 = 1;
a = 1;
w = 5;
t_end = 10;
zeta = (tau1+tau2)/(2*sqrt(tau1*tau2))

num = K*[taua 1];
den = [tau1*tau2 tau1+tau2 1];

G = tf(num,den)

t = 0:0.01:t_end;
u = a*sin(w*t);
lsim(G,u,t)