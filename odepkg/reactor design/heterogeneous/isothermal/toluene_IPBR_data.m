% Isothermal PBR

% Data
Wo = 0.0;
Wf = 10000;
Xo = 0;
yo = 1;
Po = 40;
a = 0.000098;
k = 0.00087;
Ka = 1.038;
Kc = 1.39;

yao = 0.30;
ybo = 0.45;
yco = 0.00;
ydo = 0.00;
yIo = 0.25;

Fto = 166.666;

%Preliminary Calculations
Fao = yao*Fto;
Fbo = ybo*Fto;
Fco = yco*Fto;
Fdo = ydo*Fto;
FIo = yIo*Fto;

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
