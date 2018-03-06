function [dfdt] = TICSTR_ODE(t,f)
TICSTR_data
dfdt = zeros(5,1);

%% Definição das Variáveis:
V = f(1);
Ca = f(2);
Cb = f(3);
Cc = f(4);
Cd = f(5);

%% Equações Algébricas:
ra = -(k*Ca^a*Cb^b);
rb = (b/a)*ra;
rc = -(c/a)*ra;
rd = -(d/a)*ra;
h = V/A;
Q = kv*sqrt(h);

%% Equações Diferenciais Ordinárias:
% dV/dt = IN - OUT
dfdt(1) = Qo - Q;
% dCa/dt = IN - OUT + GEN
dfdt(2)=(Cao*Qo - Ca*Q + ra*V)/V;
% dCb/dt = IN - OUT + GEN
dfdt(3)=(Cbo*Qo - Cb*Q + rb*V)/V;
% dCc/dt = IN - OUT + GEN
dfdt(4)=(Cco*Qo - Cc*Q + rc*V)/V;
% dCd/dt = IN - OUT + GEN
dfdt(5)=(Cdo*Qo - Cd*Q + rd*V)/V;

end

