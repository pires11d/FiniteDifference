% XYLOX.M
% o-Xylene oxidation PBR reactor
%
% (c)2011 Luismar M. Porto & Cíntia Soares
% Federal University of Santa Catarina - UFSC
% Florianópolis, SC, Brazil
%
% Problem definition
%
% Reaction: o-xylene + O2 (air, in excess) -- V2O5/SiO2 --> phtalic anhydride
%                  A + B                   ---------------> C
% x(1)= psi (dimensionless concentration)
% x(2)= theta (dimensionless temperature)
% x(3)= y (dimensionless pressure = P/P0
%
% xdot(1)= dpsi/dqsi
% xdot(2)= dtheta/dqsi
% xdot(3)= dy/dqsi

% ----- Part 1: Main function, cleaning up, initial conditions, call to the ODE
% solver

function xylox          % Main function

clc;                    % Clear command window
clear;                  % Clear variables and functions from memory
close all;              % Closes all the open figure windows

global L Ca0 T0 P0

% Initial conditions and integration range

ic=[1; 1; 1];                   % Initial conditions
wspan = [0 1];                  % Integration range

% Call the solver (ode45, ode15s, ...)

[w,x] = ode15s (@xylene, wspan, ic);  

% Plot

hold on
figure (1)

% Plot concentration profile

subplot (2,2,1)                 % subplot (m,n,p), m X n matrix 
plot (w*L, x(:,1)*Ca0)
title ('Concentration profile')
xlabel ('Reactor length, m')
ylabel ('Molar concentration, mol/m³')  % Cubic: Alt+Ctrl+3 

% Plot conversion

subplot (2,2,2)
plot (w*L, 1-x(:,1))
title ('Conversion')
xlabel ('Reactor length, m')
ylabel ('X, dimensionless')


% Plot temperature profile
subplot (2,2,3)
plot (w*L, x(:,2)*(T0 - 273.15))
title ('Temperature profile')
xlabel ('Reactor length, m')
ylabel ('Temperature, °C')              % Degree: Alt+Ctrl+:


% Plot pressure drop
subplot (2,2,4)
plot (w*L, x(:,3)*P0)
title ('Pressure drop')
xlabel ('Reactor length, m')
ylabel ('\Delta P, bar')

% ----- Part 2: Problem definition

function xdot = xylene (w, x) 

global L Ca0 T0 P0

%  Declaration of parameters and constants

R       = 0.082e-3;     % bar*m3/(mol*K)                    // Universal gas constant
k0      = 4.1219e8;     %                                   // Arrhenius pre-exponential factor
rhoP    = 1300;         % kg/m3                             // Catalyst particle density
epslonB = 0.38;         % dimensionless                     // Bed porosity
L       = 3;            % m                                 // Bed length
Pb0     = 0.211;        % bar                               // Inlet partial pressure for B
T0      = 352 + 273.15; % K                                 // Inlet gas temperature; Wall temperature
Ca0     = 0.02;          % mol/m3                            // Inlet concentration of A
us      = 1;            % m/s                               // Inlet superficial gas velocity
DeltaH  = 1285.409;     % kJ/mol (exothermic)               // Heat of reaction
rhoG    = 1.293;        % kg/m3                             // Inlet gas density
Cp      = 0.992;        % kJ/(kg*K)                         // gas heat capacity
Tw      = 352 + 273.15; % K                                 // Wall temperature
dt      = 2.54e-2;      % m                                 // Tube diameter
f       = 26;           % dimensionless                     // friction factor
P0      = 101325*1.50;  % Pa                                // Inlet total pressure
U       = 0.096;        % kJ/(m2*s)                         // Global heat exchange coefficient
dp      = 3e-3;         % m                                 // Particle diameter
%Mm      = 29.48;        % kg/kmol                           // Mean molecular mass

% Defined dimensionless groups

a = rhoP*(1-epslonB)*L*k0*Pb0*R*T0/us;
b = DeltaH*rhoP*L*(1-epslonB)*k0*Pb0*Ca0*R/(us*rhoG*Cp);
c = 4*U*L/(dt*us*rhoG*Cp);
d = Tw/T0;
e = f*rhoG*us^2*L/(dp*P0);

% ----- Part 3: Algebraic equations (functions of dependent variables x(i))

kappa = exp(-13636/(T0*x(2)));

% ----- Part 4: Set of dimensionless ordinary diffential equations (ODE's)

xdot(1, :)= -a*kappa*x(1)*x(2);                         % Mass (or molar) balance
xdot(2, :)= b*kappa*x(1)*x(2) + c*(d-x(2));             % Energy balance
xdot(3, :)= -e;                                         % Pressure drop (Momentum balance)

% ---- END -----