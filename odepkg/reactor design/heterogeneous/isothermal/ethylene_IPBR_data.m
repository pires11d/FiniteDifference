% Isothermal PBR

% Data
Wo = 0.0;
Wf = 45.4;
Xo = 0;
yo = 1;
Po = 10;
k = 0.0266;
Ac = 0.01414;
gc = 4.17e+8;
po = 0.413;
pc = 120;
visc = 0.0673;
Dp = 0.0208;
eb = 0.45;

Fao = 1.08;
Fbo = 0.54;
Fco = 0;
Fdo = 0;
FIo = 2.03;

MMa = 28;
MMb = 32;
MMc = 0;
MMd = 0;
MMI = 28;

a = 1;
b = 1/2;
c = 1;
d = 0;

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
