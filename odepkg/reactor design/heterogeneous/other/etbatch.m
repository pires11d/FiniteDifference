% ETBATCH.M
% Glucose Fermentation to Ethanol in a Batch Reactor
%
% (c)2009 Andrei Pavei Battisti, Jilly Anne de Souza, Veridiana Gavanski
% Federal University of Santa Catarina - UFSC
% Florian?polis, SC, Brazil
%
% Problem definition
%
% Reaction: Glucose --Saccharomyces cerevisiae--> Ethanol + More Cells
%             S     ---------------------------->   P     +     X
% x(1)= Cell Concentration (X)
% x(2)= Glucose Concentration (Substract, S)
% x(3)= Ethanol Concentration (Product, P)
% x(4)= Generation rate (rg) 
%
% xdot(1)= dX/dt
% xdot(2)= dS/dt
% xdot(3)= dP/dt


% ----- Part 1: Main function, cleaning up, initial conditions, call to the ODE
% solver

function etbatch          % Main function

clc;                    % Clear command window
clear;                  % Clear variables and functions from memory
close all;              % Closes all the open figure windows

% Initial conditions and integration range

ic=[1; 250; 0 ];               % Initial conditions
tspan = [0 12];                  % Integration range

% Call the solver (ode45, ode15s, ...)

[t,x] = ode15s (@ethanol, tspan, ic);  

% Plot

hold on
figure (1)

% Plot cells concentration profile

subplot (2,2,4)                 % subplot (m,n,p), m X n matrix 
plot (t, x(:,1),  'LineWidth',2)
axis ([0 12 0 20])
title ('Cellular concentration profile')
xlabel ('Time, h')
ylabel ('Mass concentration, g/dm³')  % Cubic: Alt+Ctrl+3 

% Plot substract concentration and Product concentration profile

subplot (2,2,2)                 % subplot (m,n,p), m X n matrix 
plot (t, x(:,2), t, x(:,3),  'LineWidth',2)
legend('Glucose', 'Ethanol');
axis ([0 12 0 300])
title ('Glucose and Ethanol concentration profile')
xlabel ('Time, h')
ylabel ('Mass concentration, g/dm³')  % Cubic: Alt+Ctrl+3 

             
% rates calculation: growth, death and maintenance
A = 1- x(:,3 )/93;
B = (A).^0.52;
C = 0.33*B.*x(:,1).*x(:,2);                          % rg
D = 1.7+ x(:,2);
RG = C./D;

RD = 0.01*x(:, 1);                                   % rd 

RSM = 0.03*x(:, 1);                                  % rsm 

% Plot rates: growth, death and cellular maintenance  
subplot (2,2,1)
plot (t, RG, 'b', t, RD, 'r', t, RSM, 'g',  'LineWidth',2)
legend('RG', 'RD','RSM');
axis ([0 12 0 2.5])
title ('Growth, death and maintenance rates')
xlabel ('Time, h')
ylabel ('Rates, g/dm³.h')

% Plot inibition by product constant
subplot (2,2,3)
KOBS = (1-(x(:,3)./93)).^0.52;
plot(t, KOBS,  'LineWidth',2)
axis ([0 12 0 1.0])
title ('Inibition by Product')
xlabel ('Time, h')
ylabel ('Kobs, dimensionless')


% ----- Part 2: Problem definition

function xdot = ethanol (t, x) 

%  Declaration of parameters and constants

Cpmax   = 93     ;      % g/dm3                             // Maximum product concentration 
n       = 0.52    ;     % dimensionless                     // Empyric constant 
Mimax   = 0.33;         % 1/h                               // Maximum growth rate 
m       = 0.03;         % gS/(gX*h)                         // Cellular maintenance
Ks      = 1.7;          % g/dm3                             // Monod constant
Kd      = 0.01;         % 1/h                               // Cellular death constant
Ysc     = 1/0.08  ;     % gX/gS                             // Yield substract/cells
Ypc     = 5.6      ;    % gP/gX                             // Yield product/cells


% ----- Part 3: Algebraic equations (functions of dependent variables x(i))


Kobs = (1-x(3)/Cpmax)^n    ;                              % Inibition by product
rg=Mimax*Kobs*x(1)*x(2)/(Ks+x(2));                        % Growth rate
rd=Kd*x(1) ;                                              % Death rate
rsm=m*x(1) ;                                              % Cellular maintenance rate

% ----- Part 4: Set of dimensionless ordinary diffential equations (ODE's)

xdot(1, :)= rg-rd      ;                              % Cellular mass balance
xdot(2, :)= Ysc*(-rg)-rsm        ;                    % Substract mass balance
xdot(3, :)= rg*Ypc      ;                             % Product balance






% ---- END -----