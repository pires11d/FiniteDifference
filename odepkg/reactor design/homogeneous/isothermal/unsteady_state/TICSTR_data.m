%% Constantes das substâncias:
MMa = 50e-3;                    %kg/mol
MMb = 50e-3;                    %kg/mol
MMc = 50e-3;                    %kg/mol
MMd = 50e-3;                    %kg/mol

%% Constantes da reação (aA + bB -> cC + dD):
a = 1;                          %
b = 2;                          %
c = 1;                          %
d = 2;                          %
Cao = 0.5;                      %mol/L
Cbo = 1.0;                      %mol/L
Cco = 0;                        %mol/L
Cdo = 0;                        %mol/L
k = 0.1;                        %s^-1*(L/mol)^(a+b-1)

%% Constantes do reator:
Vo = 0.1;                       %L
D = 5;                          %dm
H = 15;                         %dm
A = (pi*D^2)/4;                 %dm2
Vmax = A*H;                     %L
Qo = 4.0;                       %L/s
kv = 1.0;                       %L/s
rho = 1;                        %kg/L

