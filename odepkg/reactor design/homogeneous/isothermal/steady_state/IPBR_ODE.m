function dydV = IPBR_ODE(v,f)
IPBR_data

dydV = zeros(1,1);

%dXdV  = dydV(1);

x = f(1);

% Equations
Pa = Pao*(1-x)/(1+e*x);
Pb = Pao*(thetab-x)/(1+e*x);
Pc = Pao*(thetac+x)/(1+e*x);
ra = k*Pa/(1+Ka*Pa+Kb*Pb);

% Mole Balance
dydV(1) = ra/Fao;

end
