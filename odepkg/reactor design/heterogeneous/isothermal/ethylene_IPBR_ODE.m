function dydW = ethylene_IPBR_ODE(w,y)
ethylene_IPBR_data

dydW = zeros(2,1);

%dXdW  = dydW(1);
%dydW = dydW(2);

x = y(1);
y = y(2);

% Equations
G = mto/Ac;
B = (G*(1-eb)/(gc*po*Dp*eb^3))*(150*(1-eb)*visc/Dp+1.75*G)/(144*14.7);
A = 2*B/(Ac*(1-eb)*pc*Po);
ra = k*((1-x)/(1+e*x))*y;

% Mole Balance
dydW(1) = ra/Fao;

% Pressure Drop
dydW(2) = -A*(1+e*x)/(2*y);

end

