function dydz = nitrous_ODE(z,f)
nitrous_data

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

% HEAT TRANSFER FLUID 
pf = 2226-0.746*Ta;                                             % kg/m3
viscf = 4.876*0.01*0.001*exp(4680/(R*Ta));                      % Pa.s
condf = (3.745e-6*Ta-4.9249e-4)*418.4;                          % W/m.K
Cpf = (115.21-200.49e-3*Ta+126.86e-6*Ta^2)*4.184;               % J/mol.K
mf = 0.01;                                                      % kg/s
% GAS MEDIUM
pg = Po*y*10^5*MMI/(R*T);                                       % kg/m3
pg2 = Po*y2*10^5*MMI/(R*T2);                                    % kg/m3
viscg = 3.674e-7*T^0.7;                                         % Pa.s
viscg2 = 3.674e-7*T2^0.7;                                       % Pa.s
condg = 2.682e-3*(1+1.123e-3*Po*y)*T^(0.71*(1-2e-4*Po*y));      % W/m.K
condg2 = 2.682e-3*(1+1.123e-3*Po*y2)*T2^(0.71*(1-2e-4*Po*y2));  % W/m.K
Cpg = 5193;                                                     % J/kg.K
%-------------------------------------------------------------------------%

% Equations
pcz = pc*(z/zf);
eb = 0.38+0.073*(1+((D/dp-2)^2/(D/dp)^2));
Ac = pi*(D^2)/4;
G = mto/Ac;
us = G/pg;
B = (G*(1-eb)/(pg*dp*eb^3))*(150*(1-eb)*viscg/dp+1.75*G);

% Heat Transfer Coefficient
Ai = 4*pi*zf*D;
Ao = 4*pi*zf*DD;
Am = (Ai-Ao)/log(Ai/Ao);
ratio = D/DD;
Rei = (phi*dp*G)/viscg;
Ree = 4*mf/(pi*(D+DD)*viscf);
if (Ree < 2300)
    Nuo = 1.76*(3.76136-ratio);
else
    Nuo = 0.023*(Ree^(4/5))*(Cpf*viscf/condf)^0.3;
end
ho = Nuo*condf/(DD-D);
Nui = (4184/3600)*0.813*(Rei^0.9)*exp(-6*dp/D);
hi = Nui*condg/D;
U = 1/((1/hi)+(L*Ai/(condw*Am))+(Ai/(ho*Ao)));

% Reaction Rate
k = ko*exp((Ea/R)*(-1/T));
Pa = Pao*((1-x)/(1+e*x))*y;
Pb = Pao*((thetab-(b/a)*x)/(1+e*x))*y;
Pc = Pao*((thetac+(c/a)*x)/(1+e*x))*y;
Pd = Pao*((thetad+(d/a)*x)/(1+e*x))*y;
PI = Pao*((thetaI)/(1+e*x))*y;
ra = k*pcz*(1-eb)*Pb/(1+K*(Pb/Pa));

% Mole Balance
dydz(1) = 1000*ra*Ac/Fao;
% Pressure Drop
dydz(3) = -B*((1+e*x)/(Po*10^5*y))*T/To;
% Energy Balance
dydz(5) = ((-dHr*pcz*(1-eb)*ra)-(4*(U/D)*(T-Ta)))/(us*Cpg*pg);
% Catalyst Weight
dydz(7) = pcz*(1-eb)*Ac;
%-------------------------------------------------------------------------%

% Equations (PBR)
pc2 = pc*fracao_pbr;
D2 = D;
DD2 = DD;
eb2 = 0.38+0.073*(1+((D2/dp-2)^2/(D2/dp)^2));
Ac2 = pi*(D2^2)/4;
G2 = mto/Ac2;
us2 = G2/pg2;
B2 = (G2*(1-eb)/(pg2*dp*eb^3))*(150*(1-eb)*viscg2/dp+1.75*G2);

% Heat Transfer Coefficient (PBR)
Ai2 = 4*pi*zf*D2;
Ao2 = 4*pi*zf*DD2;
Am2 = (Ai2-Ao2)/log(Ai2/Ao2);
ratio2 = D2/DD2;
Rei2 = (phi*dp*G2)/viscg2;
Ree2 = 4*mf/(pi*(D2+DD2)*viscf);
if (Ree2 < 2300)
    Nuo2 = 1.76*(3.76136-ratio2);
else
    Nuo2 = 0.023*(Ree2^(4/5))*(Cpf*viscf/condf)^0.3;
end
ho2 = Nuo2*condf/(DD2-D2);
Nui2 = (4184/3600)*0.813*(Rei2^0.9)*exp(-6*dp/D2);
hi2 = Nui2*condg2/D2;
U2 = 1/((1/hi2)+(L*Ai2/(condw*Am2))+(Ai2/(ho2*Ao2)));

% Reaction Rate (PBR)
k2 = ko*exp((Ea/R)*(-1/T2));
Pa2 = Pao*((1-x2)/(1+e*x2))*y2;
Pb2 = Pao*((thetab-(b/a)*x2)/(1+e*x2))*y2;
Pc2 = Pao*((thetac+(c/a)*x2)/(1+e*x2))*y2;
Pd2 = Pao*((thetad+(d/a)*x2)/(1+e*x2))*y2;
PI2 = Pao*((thetaI)/(1+e*x2))*y2;
ra2 = k2*pc2*(1-eb2)*Pb2/(1+K*(Pb2/Pa2));

% Mole Balance (PBR)
dydz(2) = 1000*ra2*Ac2/Fao;
% Pressure Drop (PBR)
dydz(4) = -B2*((1+e*x2)/(Po*10^5*y2))*T2/To;
% Energy Balance (PBR)
dydz(6) = ((-dHr*pc2*(1-eb2)*ra2)-(4*(U2/D2)*(T2-Ta)))/(us2*Cpg*pg);
% Catalyst Weight (PBR)
dydz(8) = pc2*(1-eb2)*Ac2;
%-------------------------------------------------------------------------%

end

