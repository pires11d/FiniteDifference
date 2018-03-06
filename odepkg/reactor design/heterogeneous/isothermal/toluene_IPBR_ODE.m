function dydW = toluene_IPBR_ODE(w,f)
toluene_IPBR_data

dydW = zeros(1,1);

%dXdW  = dydW(1);

x = f(1);

% Equations
y = (1-a*w)^0.5;
Pa = Pao*(1-x)*y;
Pb = Pao*(thetab-x)*y;
Pc = Pao*(thetac+x)*y;
ra = k*Pa*Pb/(1+Ka*Pa+Kc*Pc);

% Mole Balance
dydW(1) = ra./Fao;

end

