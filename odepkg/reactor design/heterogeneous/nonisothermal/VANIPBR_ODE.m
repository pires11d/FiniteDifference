function dydz = VANIPBR_ODE(z,f)
VANIPBR_data

dydz = zeros(8,1);

%dXdz       = dydz(1);
%dXdz_pbr   = dydz(4);
%dydz       = dydz(2);
%dydz_pbr   = dydz(5);
%dTdz       = dydz(3);
%dTdz_pbr   = dydz(6);
%dWdz       = dydz(7);
%dWdz_pbr   = dydz(8);

x = f(1);
x2 = f(2);
y = f(3);
y2 = f(4);
T = f(5);
T2 = f(6);
W = f(7);
W2 = f(8);

%-------------------------------------------------------------------------%

% Equations
D = Do-(Do-Df)*((z)/zf);
Ac = pi*(D^2)/4;
G = mto/Ac;
us = G/po;
B = (G*(1-eb)/(po*Dp*eb^3))*(150*(1-eb)*visc/Dp+1.75*G);
A = 2*B/(Ac*(1-eb)*pc*Po);
k = ko*exp((Ea/R)*((1/To)-(1/T)));
Ca = Cao*((1-x)/(1+e*x))*y;
ra = k*pc*(1-eb)*Ca;

% Mole Balance
dydz(1) = ra*Ac/Fao;
% Pressure Drop
dydz(3) = -A*((1+e*x)/(2*y))*T/To;
% Energy Balance
dydz(5) = ((-dHr*pc*(1-eb)*ra)-(4*(U/D)*(T-Ta)))/(us*Cp*po);
% Catalyst Weight
dydz(7) = pc*(1-eb)*Ac;

%--------------------------------------------------------------------------%
% Equations (PBR)
D2 = (Do+Df)/2;
Ac2 = pi*(D2^2)/4;
G2 = mto/Ac2;
us2 = G2/po;
B2 = (G2*(1-eb)/(po*Dp*eb^3))*(150*(1-eb)*visc/Dp+1.75*G2);
A2 = 2*B2/(Ac2*(1-eb)*pc*Po);
k2 = ko*exp((Ea/R)*((1/To)-(1/T2)));
Ca2 = Cao*((1-x2)/(1+e*x2))*y2;
ra2 = k2*pc*(1-eb)*Ca2;

% Mole Balance (PBR)
dydz(2) = ra2*Ac2/Fao;
% Pressure Drop (PBR)
dydz(4) = -A2*((1+e*x2)/(2*y2))*T2/To;
% Energy Balance (PBR)
dydz(6) = ((-dHr*pc*(1-eb)*ra2)-(4*(U/D2)*(T2-Ta)))/(us2*Cp*po);
% Catalyst Weight (PBR)
dydz(8) = pc*(1-eb)*Ac2;
%-------------------------------------------------------------------------%

end

