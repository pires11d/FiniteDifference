% Nonisothermal Batch Reactor

% Data
ko = 1;
Ea = 10000;
Vo = 100;
R = 8.314;

a = 1;
b = 1;
c = 1;
d = 0.0;

Nao = 10;
Nbo = 10;
Nco = 0.0;
Ndo = 0.0;
NIo = 10;

Cpa = 4;
Cpb = 4;
Cpc = 4;
Cpd = 4;
CpI = 4;

dHr = -1000;
U = 1000;
A = 0.01;
To = 373.15;
Ta = 365;
Tr = 298.15;

to = 0.0;
Xf = 0.90;

%Preliminary Calculations
Nt = Nao + Nbo + Nco + Ndo + NIo;

Cao = Nao/Vo;
Cbo = Nbo/Vo;
Cco = Nco/Vo;
Cdo = Ndo/Vo;
CIo = NIo/Vo;

yao = Nao/Nt;
ybo = Nbo/Nt;
yco = Nco/Nt;
ydo = Ndo/Nt;
yIo = NIo/Nt;

thetaa = Nao/Nao;
thetab = Nbo/Nao;
thetac = Nco/Nao;
thetad = Ndo/Nao;
thetaI = NIo/Nao;

delta = (d+c-b-a)/a;
e = yao*delta;
sumCpi = Cpa*thetaa+Cpb*thetab+Cpc*thetac+Cpd*thetad+CpI*thetaI;
dCp = (d/a)*Cpd+(c/a)*Cpc-(b/a)*Cpb-Cpa;
