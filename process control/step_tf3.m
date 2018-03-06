clf
close all
clear all
clc
pkg load control

K = 1;
taua = 1;
taub = 1;
tau1 = 0.1;
tau2 = 0.1;
tau3 = 0.1;

num = K*[taua taub 1];
den = [tau1*tau2*tau3 tau1*tau2+tau1*tau3+tau2*tau3 tau1+tau2+tau3 1];
G = tf(num,den)
step(G)
