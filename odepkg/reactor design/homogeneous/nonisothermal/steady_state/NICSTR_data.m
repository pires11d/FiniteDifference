% Nonisothermal CSTR

% Data
ko = 1;
Ea = 10000;
Qo = 100;
R = 8.314;

a = 1;
b = 1;
c = 1;
d = 0.0;

Fao = 10;
Fbo = 10;
Fco = 0.0;
Fdo = 0.0;
FIo = 10;

Cpa = 4;
Cpb = 4;
Cpc = 4;
Cpd = 4;
CpI = 4;

dHr = -1000;
U = 1000;
A = 1;
To = 373.15;
Ta = 365;
Tr = 298.15;

Vo = 0.0;
Xf = 0.90;

%Preliminary Calculations
Ft = Fao + Fbo + Fco + Fdo + FIo;

Cao = Fao/Qo;
Cbo = Fbo/Qo;
Cco = Fco/Qo;
Cdo = Fdo/Qo;
CIo = FIo/Qo;

yao = Fao/Ft;
ybo = Fbo/Ft;
yco = Fco/Ft;
ydo = Fdo/Ft;
yIo = FIo/Ft;

thetaa = Fao/Fao;
thetab = Fbo/Fao;
thetac = Fco/Fao;
thetad = Fdo/Fao;
thetaI = FIo/Fao;

delta = (d+c-b-a)/a;
e = yao*delta;
sumCpi = Cpa*thetaa+Cpb*thetab+Cpc*thetac+Cpd*thetad+CpI*thetaI;
dCp = (d/a)*Cpd+(c/a)*Cpc-(b/a)*Cpb-Cpa;
