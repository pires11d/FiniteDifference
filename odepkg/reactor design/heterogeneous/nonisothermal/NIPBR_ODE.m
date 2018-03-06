function dydz = NIPBR_ODE(Z,y)

NIPBR_data

dydz = zeros(9,1);

%dCadz  = dydz(1);
%dTdz = dydz(2);
%dPdz = dydz(3);
%dPadz = dydz(4);
%dXdz = dydz(5);
%dCbdz = dydz(6);
%dPbdz = dydz(7);

Ca = y(1);
T = y(2);
P = y(3);
Pa = y(4);
X = y(5);
Cb = y(6);
Pb = y(7);
Cc = y(8);
Pc = y(9);

% Calculations
Ai = 4*pi*zf*dt^2;
Ao = 4*pi*zf*dtt^2;
Am = (Ai-Ao)/log(Ai/Ao);
e = 0.38+0.073*(1+((dt/dp-2)^2/(dt/dp)^2));
pb = (1-e)*pc;
pg = 0.001*MMo*P/(R*T);
G = pg*us;
ratio = dt/dtt;
Re = (phi*dp*G)/viscg;
f = ((1-e)/(e^3))*(1.75+150*(1-e)/Re);
Ree = 4*mf/(pi*(dt+dtt)*viscf);
if (Ree < 2300)
    Nuo = 1.76*(3.76136-ratio);
else
    Nuo = 0.023*(Ree^(4/5))*(Cpf*viscf/condf)^0.3;
end
ho = Nuo*condf/(dtt-dt);
Nui = 3.50*(Re^0.7)*exp(-4.6*dp/dt);
hi = Nui*condg/dt;
U = 1/((1/hi)+(L*Ai/(condw*Am))+(Ai/(ho*Ao)));
k = ko*exp((Ea/R)*((1/To)-(1/T)));
ra = ko*Ca/(1+Ka*Ca+Kb*Cb);

% Mole Balance of A
dydz(1) = -pb*ra/us;
% Energy Balance
dydz(2) = ((-dHr*pb*ra)-(4*(U/dt)*(T-Ta)))/(us*Cp*pg);
% Pressure Drop
dydz(3) = -f*pg*us^2/(g*dp);
% Partial Pressure of A
dydz(4) = (Pa/T)*dydz(3)-(R*T*ra*pb)/us;
% Conversion
dydz(5) = -dydz(1)/Cao;
% Mole Balance of B
dydz(6) = -pb*ra/us;
% Partial Pressure of B
dydz(7) = (Pb/T)*dydz(3)-(R*T*ra*pb)/us;
% Mole Balance of C
dydz(8) = pb*ra/us;
% Partial Pressure of C
dydz(9) = (Pc/T)*dydz(3)+(R*T*ra*pb)/us;

end

