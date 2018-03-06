function [dfdt] = TNICSTR_ODE(t,f)
TNICSTR_data
dfdt = zeros(7,1);

%% Definição das Variáveis:
V = f(1);
Ca = f(2);
Cb = f(3);
Cc = f(4);
Cd = f(5);
T = f(6);
Tc = f(7);

%% Equações Algébricas:
k = ko*exp(-Ea/(R*T));
ra = -(k*Ca^a*Cb^b);
rb = (b/a)*ra;
rc = -(c/a)*ra;
rd = -(d/a)*ra;
h = V/A;
Q = kv*sqrt(h);
Ac = pi*D*h+(pi*D^2)/4;

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
% dT/dt = IN - OUT + GEN - TRANSF
dfdt(6)=(pm*Qo*Cpm*(To-T)+dHr*ra*V-U*Ac*(T-Tc)-T*dfdt(1))/(V*Cpm*pm);
% dTc/dt = IN - OUT + TRANSF
dfdt(7)=(pf*Qco*Cpf*(Tco-Tc)+U*Ac*(T-Tc))/(Vc*Cpf*pf);

end

