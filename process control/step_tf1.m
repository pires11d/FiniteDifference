clcclf
close all
clear all
clc
pkg load control

K = 1;
tau = 1;

num = K*[1];
den = [tau 1];
G = tf(num,den)
step(G)