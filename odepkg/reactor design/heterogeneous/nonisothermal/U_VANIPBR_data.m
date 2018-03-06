% Variable Area Nonisothermal PBR

% Data
zo = 0.0;
zf = 20;
Xo = 0;
Xf = 1;
yo = 1;
Wo = 0;
Po = 500;
To = 100+273.15;
ko = 5e-8;
Ea = 10000;
R = 8.314;
Do = 0.2;
Df = 0.8;
DDo = 0.05+Do;
DDf = 0.05+Df;
L = 0.005;
po = 0.1;
pc = 1500;
phi = 1;
viscg = 1e-6;
viscf = 1e-4;
condg = 0.5;
condw = 100;
condf = 1;
mf = 10;
dp = 0.0005;
Cp = 2000;
Cpf = 5000;
dHr = -60;
Ta = 95+273.15;

Fao = 0.015;
Fbo = 0;
Fco = 0;
Fdo = 0;
FIo = 0;

MMa = 100;
MMb = 0;
MMc = 98;
MMd = 2;
MMI = 0;

a = 1;
b = 0;
c = 1;
d = 1;

%Preliminary Calculations
Fto = Fao + Fbo + Fco + Fdo + FIo;

yao = Fao/Fto;
ybo = Fbo/Fto;
yco = Fco/Fto;
ydo = Fdo/Fto;
yIo = FIo/Fto;

Pao = yao*Po;
Pbo = ybo*Po;
Pco = yco*Po;
Pdo = ydo*Po;
PIo = yIo*Po;

Cao = Pao/(R*To);
Cbo = Pbo/(R*To);
Cco = Pco/(R*To);
Cdo = Pdo/(R*To);
CIo = PIo/(R*To);

thetaa = Fao/Fao;
thetab = Fbo/Fao;
thetac = Fco/Fao;
thetad = Fdo/Fao;
thetaI = FIo/Fao;

mao = Fao*MMa*0.001;
mbo = Fbo*MMb*0.001;
mco = Fco*MMc*0.001;
mdo = Fdo*MMd*0.001;
mIo = FIo*MMI*0.001;

mto = mao + mbo + mco + mdo + mIo;
delta = (d+c-b-a)/a;
e = yao*delta;
