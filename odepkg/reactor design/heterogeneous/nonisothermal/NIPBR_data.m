% Nonisothermal PBR with pressure drop

% Data
R = 8.314;
To = 100+273.15;
Po = 1000;
Xo = 0;
yao = 1/3;
ybo = 2/3;
yco = 0;
Pao = yao*Po;
Pbo = ybo*Po;
Pco = yco*Po;
Cao = Pao/(R*To);
Cbo = Pbo/(R*To);
Cco = Pco/(R*To);
zf = 2;

a = 1;
b = 1;
c = 1;
d = 0.0;

ko = 1e-3;
Ea = 10000;
Ka = 0.05;
Kb = 0.1;
Cp = 1000;
Cpf = 5000;
dHr = -20000;
Ta = 70+273.15;

phi = 1;
pc = 1300;
us = 0.5;
dt = 0.100;
dtt = 0.120;
dp = 0.005;
MMo = 100;
g = 9.81;
mf = 1;
L = 0.005;
viscg = 0.0001;
viscf = 0.001;
condg = 0.02;
condw = 100;
condf = 0.58;

