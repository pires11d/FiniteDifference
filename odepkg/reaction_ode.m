function dCdx = reaction_ode(x,C)
dCdx = zeros(6,1);

reaction_data

%% Definitions
Ca = C(1);
dCa = C(2);
Cb = C(3);
dCb = C(4);
Cr = C(5);
dCr = C(6);

r1 = k1*Ca*Cb;
r2 = k2*Ca*Cr;

% ODEs
dCdx(1) = dCa;
dCdx(2) = (1/Da)*(r1+r2);
dCdx(3) = dCb;
dCdx(4) = (1/Db)*r1;
dCdx(5) = dCr;
dCdx(6) = (1/Dr)*(r2-r1);
end
