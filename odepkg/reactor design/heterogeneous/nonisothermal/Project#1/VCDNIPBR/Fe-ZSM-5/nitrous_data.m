% Lucas Joshua Pires
% 12100929
% Departamento de Eng. Química e Eng. de Alimentos - EQA
% Universidade Federal de Santa Catarina - UFSC
%-------------------------------------------------------------------------%
%------------ NONISOTHERMAL HETEROGENEOUS PACKED-BED REACTOR -------------%
%--------- with catalyst density change in the axial direction -----------%
%-------------------------------------------------------------------------%
%----------------------------- General Data ------------------------------%
%-------------------------------------------------------------------------%
% INITIAL/FINAL CONDITIONS:
zo = 0.0;           % m
zf = 0.5;           % m
Xo = 0;             % dimensionless
Xf = 1;             % dimensionless
yo = 1;             % dimensionless
Wo = 0;             % kg
Po = 1;             % bar
To = 640+273.15;    % K
Ta = To-10;         % K
% STOICHIOMETRY:
a = 1;              % N2O
b = 1;              % NO
c = 1;              % N2
d = 0;              % O2
dd = 1;             % NO2
% MOLECULAR WEIGHTS:
MMa = 44e-3;        % kg/mol N2O
MMb = 0;            % kg/mol NO
MMc = 28e-3;        % kg/mol N2
MMd = 32e-3;        % kg/mol O2
MMdd = 46e-3;       % kg/mol NO2
MMI = 4.02e-3;      % kg/mol He
% FEED CONDITIONS:
Qo = 1e-6;          % m3/s
Pao = 0.030;        % bar
Pbo = 0.060;        % bar
Pco = 0.000;        % bar
Pdo = 0.000;        % bar
Pddo = 0.000;       % bar
PIo=Po-Pao-Pbo-Pdo; % bar
% REACTION DATA:
R = 8.314;          % J/mol.K
dHr = -81500;       % J/mol
%---Co-ZSM-5
%     ko = 1.7e3;     % mol/s.bar.kg                
%     Ea = 104000;    % J/mol
%     K = 4.5;        % dimensionless
%---Cu-ZSM-5
%     ko = 2.04e3;    % mol/s.bar.kg                
%     Ea = 136000;    % J/mol
%     K1 = 413;       % /bar
%     K2 = 363;       % /bar
%---Fe-ZSM-5
    ko = 98e3;      % mol/s.bar.kg                
    Ea = 168000;    % J/mol
    K = 14;         % dimensionless
% REACTOR VESSEL:
D = 0.010;          % m
DD = D*1.5;         % m
fracao_pbr = 2/4;   % dimensionless
L = 0.002;          % m
condw = 230;        % W/m.K
% CATALYST PELLETS:
phi = 1;            % sphericity
dp = 0.000150;      % m
pc = 3220;          % kg/m3

%-------------------------------------------------------------------------%
%------------------------ Preliminary Calculations -----------------------%
%-------------------------------------------------------------------------%
% Initial Molar Fraction
yao = Pao/Po;       
ybo = Pbo/Po;
yco = Pco/Po;
ydo = Pdo/Po;
yddo = Pddo/Po;
yIo = PIo/Po;
% Initial Pressure Ratio
thetaa = Pao/Pao;
thetab = Pbo/Pao;
thetac = Pco/Pao;
thetad = Pdo/Pao;
thetadd = Pddo/Pao;
thetaI = PIo/Pao;
% Initial Concentration
Cao = 10^5*Pao/(R*To);      % mol/m3
Cbo = 10^5*Pbo/(R*To);      % mol/m3
Cco = 10^5*Pco/(R*To);      % mol/m3
Cdo = 10^5*Pdo/(R*To);      % mol/m3
Cddo = 10^5*Pddo/(R*To);    % mol/m3
CIo = 10^5*PIo/(R*To);      % mol/m3
% Initial Molar Flow rate
Fao = Cao*Qo;               % mol/s
Fbo = Cbo*Qo;               % mol/s
Fco = Cco*Qo;               % mol/s
Fdo = Cdo*Qo;               % mol/s
Fddo = Cddo*Qo;             % mol/s
FIo = CIo*Qo;               % mol/s
% Initial Mass Flow rate
mao = Fao*MMa;              % kg/s
mbo = Fbo*MMb;              % kg/s
mco = Fco*MMc;              % kg/s
mdo = Fdo*MMd;              % kg/s
mddo = Fddo*MMdd;           % kg/s
mIo = FIo*MMI;              % kg/s
% Initial Total Molar Flow rate
Fto = Fao + Fbo + Fco + Fdo + Fddo + FIo;
% Initial Total Mass Flow rate
mto = mao + mbo + mco + mdo + mddo + mIo;
% Expansion Coefficient 
delta = (dd+c-b-a)/a;
e = yao*delta;
