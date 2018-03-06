% Nonisothermal PFR

% Data
Vo = 0.0;
To = 1035;
Po = 162;
Xo = 0;
Xf = 0.99;

dHr = 80770;
U = 110;
A = 150;
Ta = 1150;
Tr = 298.15;
R = 8.314;

Fao = 38.3;
Fbo = 0;
Fco = 0;
Fdo = 0;
FIo = 0;

a = 1;
b = 0;
c = 1;
d = 1;

%Preliminary Calculations
Ft = Fao + Fbo + Fco + Fdo + FIo;

yao = Fao/Ft;
ybo = Fbo/Ft;
yco = Fco/Ft;
ydo = Fdo/Ft;
yIo = FIo/Ft;

Pao = yao*Po;
Pbo = ybo*Po;
Pco = yco*Po;
Pdo = ydo*Po;
PIo = yIo*Po;

Cao = 1000*Pao/(R*To);
Cbo = 1000*Pbo/(R*To);
Cco = 1000*Pco/(R*To);
Cdo = 1000*Pdo/(R*To);
CIo = 1000*PIo/(R*To);

thetaa = Fao/Fao;
thetab = Fbo/Fao;
thetac = Fco/Fao;
thetad = Fdo/Fao;
thetaI = FIo/Fao;

delta = (d+c-b-a)/a;
e = yao*delta;
