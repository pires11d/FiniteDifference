clf
close all
clear all
clc

digits(2)
C0 = zeros(1,4);
C = fsolve('CSTR_f',C0);
clc
Ca = C(1)
Cb = C(2)
Cc = C(3)
Cd = C(4)







