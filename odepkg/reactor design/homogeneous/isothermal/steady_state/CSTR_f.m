function f = CSTR_f(C)
k1 = 1.5;
k2 = 0.1;
k3 = 0.1;
k4 = 0.5;
Q = 50;
V = 100;
Cao = 1;
Cbo = 0;
Cco = 0;
Cdo = 0;

Co = [Cao;Cbo;Cco;Cco];

f(1)=-C(1)+Co(1)+V*((-k1*C(1))-k2*C(1)^1.5+k3*C(3)^2)/Q;
f(2)=-C(2)+Co(2)+V*(2*k1*C(1)-k4*C(2)^2)/Q;
f(3)=-C(3)+Co(3)+V*(k2*C(1)^1.5-k3*C(3)^2+k4*C(2)^2)/Q;
f(4)=-C(4)+Co(4)+V*(k4*C(2)^2)/Q;
end

