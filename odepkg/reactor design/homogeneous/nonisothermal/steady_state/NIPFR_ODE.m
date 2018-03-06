function dydX = NIPFR_ODE(X,y)
NIPFR_data

dydX = zeros(2,1);

%dVdX  = dydX(1);
%dTdX = dydX(2);

V = y(1);
T = y(2);

% Reaction Rate
k = ko*exp((Ea/R)*((1/To)-(1/T)));
Ca = Cao*((1-X)/(1+e*X))*(To/T);
ra = k*Ca;

% Mole Balance
dydX(1) = Fao/ra;

% Energy Balance
dydX(2) = (U*A*(Ta-T)+ra*(-dHr+dCp*(T-To)))/(ra*(sumCpi+X*dCp));

end

