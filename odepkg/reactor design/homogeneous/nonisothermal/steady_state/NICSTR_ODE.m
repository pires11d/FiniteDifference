function y = NICSTR_ODE(V,T)
NICSTR_data

y = zeros(2,1);

V = y(1);
T = y(2);

% Reaction Rate
k = ko*exp((Ea/R)*((1/To)-(1/T)));
Ca = Cao*((1-X)/(1+e*X))*(To/T);
ra = k*Ca;

% Mole Balance
V = Fao*X/ra;

% Energy Balance
T = (U*A*(Ta-T)+ra*(-dHr+dCp*(T-To)))/(ra*(sumCpi+X*dCp));

end

