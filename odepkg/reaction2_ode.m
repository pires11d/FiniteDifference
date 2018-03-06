function dCdx = reaction2_ode(x,C)
dCdx = zeros(6,1);

reaction2_data

%% Definitions
Ca = C(1);
dCa = C(2);
Cb = C(3);
dCb = C(4);
Cr = C(5);
dCr = C(6);

%r1 = k1*Ca*Cb;
%r2 = k2*Ca*Cr;

% ODEs
dCdx(1) = dCa;
dCdx(2) = (Haa^2)*Ca*Cb*(S+1)/S;
dCdx(3) = dCb;
dCdx(4) = (Hab^2)*Ca*Cb;
dCdx(5) = dCr;
dCdx(6) = (Har^2)*Ca*Cr*(1-S);
end
