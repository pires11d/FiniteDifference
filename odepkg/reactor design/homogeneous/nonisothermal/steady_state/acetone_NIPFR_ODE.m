function dydX = NIPFR_ODE(X,y)
acetone_NIPFR_data

dydX = zeros(2,1);

%dVdX  = dydX(1);
%dTdX = dydX(2);

V = y(1);
T = y(2);

% Specific Heats
Cpa = 26.63+0.183*T-45.86*(10^-6)*T^2;
Cpb = 0;
Cpc = 20.04+0.0945*T-30.95*(10^-6)*T^2;
Cpd = 13.39+0.077*T-18.71*(10^-6)*T^2;
CpI = 0;

sumCpi = Cpa*thetaa+Cpb*thetab+Cpc*thetac+Cpd*thetad+CpI*thetaI;
dCp = (d/a)*Cpd+(c/a)*Cpc-(b/a)*Cpb-Cpa;

% Reaction Rate
k = 3.58*exp(34.222*(1/1035-1/T));
Ca = Cao*((1-X)/(1+e*X))*(To/T);
ra = k*Ca;

% Mole Balance
dydX(1) = Fao/ra;

% Energy Balance
dydX(2) = (U*A*(Ta-T)+ra*(-dHr+dCp*(T-To)))/(ra*(sumCpi+X*dCp));

end

