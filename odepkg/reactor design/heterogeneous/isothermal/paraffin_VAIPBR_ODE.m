function dydz = paraffin_VAIPBR_ODE(z,f)
paraffin_VAIPBR_data

dydz = zeros(4,1);

%dXdz = dydz(1);
%dydz = dydz(2);
%dXdz_pbr  = dydz(3);
%dydz_pbr = dydz(4);

x = f(1);
y = f(2);
xi = f(3);
yi = f(4);

% Equations
D = Do-(Do-Df)*((zf-z)/zf);
Ac = pi*(D^2)/4;
G = mto/Ac;
B = (G*(1-eb)/(po*Dp*eb^3))*(150*(1-eb)*visc/Dp+1.75*G)*(0.01);
A = 2*B/(Ac*(1-eb)*pc*Po);
Ca = Cao*((1-x)/(1+e*x))*y;
ra = k*pc*(1-eb)*Ca;

% Mole Balance
dydz(1) = ra*Ac/Fao;

% Pressure Drop
dydz(2) = -A*(1+e*x)/(2*y);

%--------------------------------------------------------------------------%
% Equations (PBR)
Di = (Do+Df)/2;
Aci = pi*(Di^2)/4;
Gi = mto/Aci;
Bi = (Gi*(1-eb)/(po*Dp*eb^3))*(150*(1-eb)*visc/Dp+1.75*Gi)*(0.01);
Ai = 2*Bi/(Aci*(1-eb)*pc*Po);
Cai = Cao*((1-xi)/(1+e*xi))*yi;
rai = k*pc*(1-eb)*Cai;

% Mole Balance (PBR)
dydz(3) = rai*Aci/Fao;

% Pressure Drop (PBR)
dydz(4) = -Ai*(1+e*xi)/(2*yi);

end

