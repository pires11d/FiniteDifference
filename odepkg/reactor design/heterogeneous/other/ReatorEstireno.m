% Universidade Federal de Santa Catarina
% Centro Tecnol�gico
% Departamento de Engenharia Qu�mica e Engenharia de Alimentos
% C�lculo de Reatores II
% Data: 02/12/15
% Professor Dr. Luismar Marques Porto
% Alunos: Bianca Machado, Matheus Tanaka, Joice Sapatieri
% Produ��o de Estireno atrav�s da desidrogena��o de Etilbenzeno
 
%%
 
function ReatorEstireno
close;
close all;
clear;
clearvars -global;
clc;
 
global T ConcA ConcB ConcC  L Ka Kb Kc R Pi DeltaH Cp U M rhob Eb dt dp mi Tw FA conve
 
% Vari�veis
 
T = 850;                                  % Temperatura absoluta [K] OK
ConcA = 10;                               % Concentra��o total de entrada [mol/m^3] DEFINIR
 
% Base de Dados
 
Ka = 2.213e-4*exp(65830*(R*T)^-1);        % Constante de adsor��o do Etilbenzeno
Kb = 8.839e-14*exp(209396*(R*T)^-1);      % Constante de adsor��o do Estireno 
Kc = 3.625e-7*exp(103151*(R*T)^-1);       % Constante de adsor��o do Hidrog�nio


Tw = 800;                                 % Temperatura da parede [K] DEFINIR
Pi = 0.4e6;                               % Press�o inicial [Pa] DEFINIR
R = 8.31446;                              % Constante universal dos gases [J/(mol.K)] 
L = 3;                                  % Comprimento do tubo [m] DEFINIR
DeltaH = 117600;                          % Calor de rea��o [J/mol] OK
Cp = 2656;                                % Capacidade calor�fica m�dia do Etilbenzeno [J/(Kg.K)] OK
U = 210;                                  % Coeficiente global de transfer�ncia de calor [J/(m^2.s)] VERIFICAR COMO CALCULAR
dt = 0.003 ;                                % Di�metro do tubo [m] DEFINIR
dp = 0.0001;                              % Di�metro da part�cula de catalisador [m] DEFINIR
mi = 6.69*10^-4;                          % Viscosidade do fluido [Pa/s)ENCONTRAR
FA = 1000;                                % Produ��o de Estireno [Kg/h] DEFINIR
M = 5.87*10^-3;                           % Massa molar m�dia [kg/mol] ENCONTRAR
rhob = 2146.3;                            % Massa espec�fica do catalisador [Kg/m^3] Fe2O3
Eb = 0.5;                                 % Porosidade do bulk DEFINIR
ConcB = 0.000001;
ConcC = 0.000001;
%-----------------------------------------%
      
z0=[T; Pi; ConcA; ConcB; ConcC];          % Condi��es iniciais   
zspan = [0 L];                            % Intervalo de integra��o       
 
sol = ode45 (@reator, zspan, z0); 
 
sol;
 

% Vetor solu��o
 
u1 = sol.y(1,:);                          % Temperatura
u2 = sol.y(2,:);                          % Press�o
u3 = sol.y(3,:);                          % Concentra��o de Etilbenzeno [mol/m^3]
u4 = sol.y(4,:);                          % Concentra��o de Estireno [mol/m^3]
u5 = sol.y(5,:);                          % Concentra��o de Hidrog�nio [mol/m^3]
z = sol.x(:);
 
 
% Plotagem dos gr�ficos


  
figure (1)

plot(z,u3,'r')    
title('Concentra��o Etilbenzeno')
xlabel('Dist�ncia (m)')
ylabel('[EB] (mol/m^3)')

figure(2)
plot(z,u4,'g')    
title('Concentra��o do Estireno')
xlabel('Dist�ncia (m)')
ylabel('[Est] (mol/m^3)')


figure(3)
plot(z,u1)  
title('Temperatura')
xlabel('Dist�ncia (m)')
ylabel('T (K)')


figure(4)
plot(z,u2)  
title('Perda de press�o')
xlabel('Dist�ncia (m)')
ylabel('P (Pa)')

figure(5)
plot(z,conve)  
title('Convers�o')
xlabel('Dist�ncia (m)')
ylabel('x')


 
 
function xD = reator (z,x)
global ConcA conve M PB k PC PA Eb DeltaH Cp U dt dp mi G Ka Kb Kc K  R rhob fp Tw FA Cb
 
% Base de Dados
Ka = 2.213e-4*exp(65830/(R*x(1)));                  % Constante de adsor��o do Etilbenzeno
Kb = 8.839e-14*exp(209396/(R*x(1)));                % Constante de adsor��o do Estireno 
Kc = 3.625e-7*exp(103151/(R*x(1)));                 % Constante de adsor��o do Hidrog�nio
K = exp(-18.7+0.021*x(1));                          
k = 5.37e-1*exp(-93400/(R.*x(1)*3600*101325));      % Constante cin�tica da rea��o [mol/(J.kg)^-1]
DeltaH = 117600;                                    % Calor de rea��o [J/mol]
PA = x(3)*R*x(1);                                   % Press�o parcial de Etilbenzeno [Pa]
PB = x(4)*R*x(1);                                   % Press�o parcial de Estireno [Pa]
PC = x(5)*R*x(1);                                   % Press�o parcial de Hidrog�nio [Pa]
rhog = x(2)*M/(R*x(1));                             % Massa espec�fica do fluido reagente [kg/m^3]
rhob = 2146.3;                                      % Massa espec�fica do leito
Eps = 1;                                            % Coeficiente de expans�o/contra��o
conve = (ConcA-x(3))/ConcA;                         % Convers�o
G = FA/conve;                                       % Fluxo m�ssico de Etilbenzeno
us = 1.5*(1+Eps*conve);                             % Velocidade superficial do fluido [m/s]
fp = (1-Eb)^2*150/(Eb^3*dp.*G/mi);                  % Fator de atrito
Denom = (1 + Ka*PA + Kb*PB + Kc*PC);                % Denominador 
rate = (k*((Ka*PA-(1/K)*Kb*Kc*PB*PC)))/Denom^2;     % Velocidade de rea��o




xD(1,:) = -1*DeltaH*rhob*(1-Eb)*rate/(us*rhog*Cp) + 4*U.*(Tw-x(1))/(dt*Cp*us*rhog);
xD(2,:) = -fp*rhog*us^2/dp;
xD(3,:) = -1*rhob*(1-Eb)*rate/us;
xD(4,:) = rhob*(1-Eb)*rate/us;
xD(5,:) = rhob*(1-Eb)*rate/us;

%-----------------------------------------%
%-----------------------------------------%
