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
ko = 0.0005;
Ea = 10000;
R = 8.314;
Do = 0.2;
Df = 0.8;
po = 1;
pc = 1500;
visc = 1e-6;
Dp = 0.001;
eb = 0.53;
Cp = 2000;
dHr = -60;
U = 20;
Ta = 95+273.15;

Fao = 0.05;
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
