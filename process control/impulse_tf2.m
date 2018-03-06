clf
close all
clear all
clc
pkg load control

K = 1;
taua = 1;
tau1 = 2;
tau2 = 4;

num = K*[taua 1];
den = [tau1*tau2 tau1+tau2 1];
G = tf(num,den)
impulse(G)
