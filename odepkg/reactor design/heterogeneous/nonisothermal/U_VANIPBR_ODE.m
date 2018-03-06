function dydz = U_VANIPBR_ODE(z,f)
U_VANIPBR_data

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
D = Do-(Do-Df)*(z/zf);
DD = DDo-(DDo-DDf)*(z/zf);
height = (Df-Do)/2;
diagonal = (zf^2+height^2)^0.5;
eb = 0.38+0.073*(1+((D/dp-2)^2/(D/dp)^2));
Ac = pi*(D^2)/4;
G = mto/Ac;
us = G/po;
B = (G*(1-eb)/(po*dp*eb^3))*(150*(1-eb)*viscg/dp+1.75*G);
A = 2*B/(Ac*(1-eb)*pc*Po);

% Heat Transfer Coefficient
Ai = 4*pi*diagonal*(Do+Df)/2;
Ao = 4*pi*diagonal*(DDo+DDf)/2;
Am = (Ai-Ao)/log(Ai/Ao);
ratio = D/DD;
Rei = (phi*dp*G)/viscg;
%f = ((1-eb)/(eb^3))*(1.75+150*(1-eb)/Rei);
Ree = 4*mf/(pi*(D+DD)*viscf);
if (Ree < 2300)
    Nuo = 1.76*(3.76136-ratio);
else
    Nuo = 0.023*(Ree^(4/5))*(Cpf*viscf/condf)^0.3;
end
ho = Nuo*condf/(DD-D);
Nui = 3.50*(Rei^0.7)*exp(-4.6*dp/D);
hi = Nui*condg/D;
U = 1/((1/hi)+(L*Ai/(condw*Am))+(Ai/(ho*Ao)));

% Reaction Rate
k = ko*exp((Ea/R)*((1/To)-(1/T)));
Pa = Pao*((1-x)/(1+e*x))*y;
ra = k*pc*(1-eb)*Pa;

% Mole Balance
dydz(1) = ra*Ac/Fao;
% Pressure Drop
dydz(3) = -A*((1+e*x)/(2*y))*T/To;
% Energy Balance
dydz(5) = ((-dHr*pc*(1-eb)*ra)-(4*(U/D)*(T-Ta)))/(us*Cp*po);
% Catalyst Weight
dydz(7) = pc*(1-eb)*Ac;
%-------------------------------------------------------------------------%

% Equations (PBR)
D2 = (Do+Df)/2;
DD2 = (DDo+DDf)/2;
eb2 = 0.38+0.073*(1+((D2/dp-2)^2/(D2/dp)^2));
Ac2 = pi*(D2^2)/4;
G2 = mto/Ac2;
us2 = G2/po;
B2 = (G2*(1-eb)/(po*dp*eb^3))*(150*(1-eb)*viscg/dp+1.75*G2);
A2 = 2*B2/(Ac2*(1-eb)*pc*Po);

% Heat Transfer Coefficient (PBR)
Ai2 = 4*pi*zf*D2;
Ao2 = 4*pi*zf*DD2;
Am2 = (Ai2-Ao2)/log(Ai2/Ao2);
ratio2 = D2/DD2;
Rei2 = (phi*dp*G2)/viscg;
%f2 = ((1-eb2)/(eb2^3))*(1.75+150*(1-eb2)/Rei2);
Ree2 = 4*mf/(pi*(D2+DD2)*viscf);
if (Ree2 < 2300)
    Nuo2 = 1.76*(3.76136-ratio2);
else
    Nuo2 = 0.023*(Ree2^(4/5))*(Cpf*viscf/condf)^0.3;
end
ho2 = Nuo2*condf/(DD2-D2);
Nui2 = 3.50*(Rei2^0.7)*exp(-4.6*dp/D2);
hi2 = Nui2*condg/D2;
U2 = 1/((1/hi2)+(L*Ai2/(condw*Am2))+(Ai2/(ho2*Ao2)));

% Reaction Rate (PBR)
k2 = ko*exp((Ea/R)*((1/To)-(1/T2)));
Pa2 = Pao*((1-x2)/(1+e*x2))*y2;
ra2 = k2*pc*(1-eb2)*Pa2;

% Mole Balance (PBR)
dydz(2) = ra2*Ac2/Fao;
% Pressure Drop (PBR)
dydz(4) = -A2*((1+e*x2)/(2*y2))*T2/To;
% Energy Balance (PBR)
dydz(6) = ((-dHr*pc*(1-eb2)*ra2)-(4*(U2/D2)*(T2-Ta)))/(us2*Cp*po);
% Catalyst Weight (PBR)
dydz(8) = pc*(1-eb2)*Ac2;
%-------------------------------------------------------------------------%

end

