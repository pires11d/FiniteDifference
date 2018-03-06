clf
close all
clear all
clc
pkg load control

K = 1;
tau = 1;
a = 1;
w = 5;
t_end = 10;

num = K*[1];
den = [tau 1];

G = tf(num,den)

t = 0:0.01:t_end;
u = a*sin(w*t);
lsim(G,u,t)
