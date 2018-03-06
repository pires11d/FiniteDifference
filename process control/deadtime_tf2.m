clf
close all
clear all
clc
pkg load control

K = 1;
taua = 0;
tau1 = 1;
tau2 = 1;
theta = 1;
zeta = (tau1+tau2)/(2*sqrt(tau1*tau2))

s = tf('s');
num = K*[taua 1];
den = [tau1*tau2 tau1+tau2 1];

G = tf(num,den)
step(G*exp(-theta*s))