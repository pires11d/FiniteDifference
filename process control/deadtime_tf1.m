clf
close all
clear all
clc
pkg load control

K = 1;
tau = 1;
theta = 1;

s = tf('s');
num = K*[1];
den = [tau 1];

G = tf(num,den)
step(G*exp(-theta*s))
