% Variable Area Isothermal PBR

% Data
zo = 0.0;
zf = 20;
Xo = 0;
Xf = 1;
yo = 1;
Po = 2000;
k = 0.02;
Do = 8;
Df = 5;
po = 0.032;
pc = 2.6;
visc = 1.5e-6;
Dp = 0.02;
eb = 0.4;

Fao = 0.44;
Fbo = 0;
Fco = 0;
Fdo = 0;
FIo = 0;

Cao = 0.32;
Cbo = 0;
Cco = 0;
Cdo = 0;
CIo = 0;

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

thetaa = Fao/Fao;
thetab = Fbo/Fao;
thetac = Fco/Fao;
thetad = Fdo/Fao;
thetaI = FIo/Fao;

mao = Fao*MMa;
mbo = Fbo*MMb;
mco = Fco*MMc;
mdo = Fdo*MMd;
mIo = FIo*MMI;

mto = mao + mbo + mco + mdo + mIo;
delta = (d+c-b-a)/a;
e = yao*delta;
