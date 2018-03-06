% XYLOX.M
% o-Xylene oxidation PBR reactor
%
% (c)2008 Luismar M. Porto
% Federal University of Santa Catarina - UFSC
% Florianopolis, SC, Brazil
% Edited in May/09 by Brunno Bagnariolli
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


% ----- Part 1: 
%       Main function, cleaning up, initial conditions, 
%       call to the ODE solver

function xyloxU          % Main function

clc;                    % Clear command window
clear all;              % Clear variables and functions from memory
close all;              % Closes all the open figure windows

global L Ca0 T0 P0

% Initial conditions and integration range

ic=[1; 1; 1];            % Initial conditions
wspan = [0 1];           % Integration range

% Call the solver (ode45, ode15s, ...)

[w,x] = ode15s (@xylene, wspan, ic);  

% Plot
hold on
figure (1)

% Plot concentration profile

subplot (2,2,1)
plot (w*L, x(:,1)*Ca0)
title ('Concentration profile')
xlabel ('Reactor length, m')
ylabel ('Molar concentration, mol/m³')

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
ylabel ('Temperature, ºC')


% Plot pressure drop
subplot (2,2,4)
plot (w*L, x(:,3)*P0)
title ('Pressure drop')
xlabel ('Reactor length, m')
ylabel ('\Delta P, bar')


% ----- Part 2: 
%       Problem definition


function xdot = xylene (w, x) 

global L Ca0 T0 P0 Tw dp dt rhoG CpG us T f

%  Declaration of parameters and constants


R       = 0.082e-3;           % bar*m3/(mol*K)                    // Universal gas constant
k0      = 4.1219e8;           %                                   // Arrhenius pre-exponential factor
rhoP    = 1300;               % kg/m3                             // Catalyst particle density
L       = 3;                  % m                                 // Bed length
Pb0     = 0.211;              % bar                               // Inlet partial pressure for B
T0      = 352 + 273.15;       % K                                 // Inlet gas temperature
Ca0     = 0.02;               % mol/m3                            // Inlet concentration of A
Cb0     = Pb0/(R*T0);         % mol/m3                            // Inlet concentration of B
us      = 1;                  % m/s                               // Inlet superficial gas velocity
DeltaH  = 1285.409;           % kJ/mol (exothermic)               // Heat of reaction
Tw      = 352 + 273.15;       % K                                 // Wall temperature
dt      = 2.54e-2;            % m                                 // Tube diameter
P0      = 101325*1.50;        % Pa                                // Inlet total pressure
dp      = 3e-3;               % m                                 // Particle diameter
EoverR  = -13636;             % K                                 // Energy of activation over universal gas constant


% ----- Part 3: Algebraic equations (functions of T)


T       = x(2)*T0;                  % K                           // Gas temperature
k       = k0*exp(EoverR/T);         %                             // Arrhenius equation
U       = globalheat;               % kJ/(m2*s)                   // Global heat exchange coefficient


% Defined dimensionless groups


alfa    = -rhoP * L * k * Cb0 * R^2 * T0^2 / us;
beta    = DeltaH * rhoP * L * k * Cb0 * R^2 * T0 * Ca0 / (us * rhoG * CpG);
gama    = 4 * U * L / (dt * us * rhoG * CpG);
delta   = -f * rhoG * us^2 * L / (dp*P0);


% ----- Part 4: Set of dimensionless ordinary diffential equations (ODE's)


xdot(1,:) = alfa * x(1) * x(2);
xdot(2,:) = beta * x(1) * x(2) - gama * (x(2) - Tw/T0);
xdot(3,:) = delta;


% ----- Complementary function:
%       Calculus of global heat exchange coefficient
%       Added in May/09 by Brunno Bagnariolli

function U = globalheat

global  Tw dp dt rhoG CpG us T f

% ----- Mixture definition

Tms     =   320;                                                % K                     // Temperature of molten salt
mxy     =   0.01;                                               % dimensionless         // Mass percent of o-xylene in the gas mixture
mair    =   1-mxy;                                              % dimensionless         // Mass percent of air in the gas mixture

  
% ----- PIPE

% Pipe Size
% Taken from Wikipedia for NPS 1
% Data for tther schedules and NPS can be found at the same website
% Select the desired pipe schedule and uncomment (Ctrl+T) the line

dte=33.4e-3;                                                    % m                 // Outside diameter for NPS 1

sch     =   [dte-1.651e-3;1.651e-3];                             % Schedule 5
% sch     =   [dte-2.769e-3;2.769e-3];                           % Schedule 10
% sch     =   [dte-2.896e-3;2.896e-3];                           % Schedule 20
% sch     =   [dte-3.378e-3;3.378e-3];                           % Schedule 40

dti     =   sch(1);                                             % m                 // Inside diameter
d       =   sch(2);                                             % m                 // Wall thickness

% Thermal conductivity for pipe wall. 
% Taken from book Fundamentals of Heat and Mass Transfer - Incropera and DeWitt
% Correlation for other materials can be found at the same book

lambdap =   1.25073443246*Tw^0.433569017045;                    % W/m*K                 // 304 Stainless Steel pipe (100K to 1000K)
% lambdap =   74.59 - 0.04445*Tw;                               % W/m*K                 // Carbon Steel (400K to 800K)


% ----- MOLTEN SALT

% Taken from various articles indexed in Tohoku Molten Salt Database
% http://ras.material.tohoku.ac.jp/~molten/molten_ps_query1.php
% Correlation for other materials can be found at the same website
% Select the desired molten salt and uncomment (Ctrl+T) the line
% The properties in matrix are:
% [ Thermal conductivity (W/m*K) ; Specific Mass (kg/m³) ; Viscosity(Pa*s) ; Cp(kJ/kg*J) ]

% salt    =   [ 0.18 + 0.000404*Tms ; ...
%             (2.3063 - 0.0007235*Tms)*1E3 ; ...
%             (29.7085 - 0.0711208*Tms + 4.47023E-5*Tms^2)*1E-3 ; ...
%             1450.70 ];                                                                % KNO3

salt    =   [ 0.28 + .0005*Tms ; ...
            (2.116 - 0.000729*Tms)*1E3 ; ...
            (0.09139*exp(4037/(8.314*Tms)))*1E-3 ; ...
            1537.99 ];                                                                  % NaNO3
    

lambdams     =   salt(1);                                       % W/m*K                 // Thermal conductivity of molten salt
rhoms        =   salt(2);                                       % kg/m³                 // Specific Mass of molten salt
mims         =   salt(3);                                       % Pa*s                  // Viscosity of molten salt
cpms         =   salt(4);                                       % kJ/kg*J               // Cp of molten salt


% ----- GAS

% Thermal conductivity
% Taken from the DIPPR database

lambdaair    =   1.5207E-11*T^3 - 4.8574E-08*T^2 ...
                + 1.0184E-04*T - 3.9333E-04;                    % W/m*K                 //  Thermal conductivity of air
lambdaxylene =   1.9989E-1-2.2990E-4*T;                         % W/m*K                 //  Thermal conductivity of o-xylene

lambdag      =   mair*lambdaair+mxy*lambdaxylene;               % W/m*K                 //  Thermal conductivity of the gas


% Viscosity

miair        =   (4.38E-2*(1.01E-3*T-0.093)^(5/9))*1E-3;         % Pa*s                  // Viscosity of Air (250K to 2500K) - From H. Y. Lo, "Viscosity of Gaseous Air at Moderate and High Pressures", J. Chem. Eng. Data, 1966, 11 (4), 540-544
mixylene     =   (3.8077E-06*T^3.1520E-01)/(1+7.7431E+02/T);     % Pa*s                  // Viscosity of o-xylene - From the DIPPR database

mig          =   mair*miair+mxy*mixylene;                        % Pa*s                  // Viscosity of the gas

% Specific Mass
% Taken from the DIPPR database

rhoG         =     361.77819*T^-1.00336;                         % kg/m³                // Specific mass of Air

% Heat Capacity
% Taken from the DIPPR database

Cpair     =     (1.9327E-10*T^4 ...
                - 7.9999E-07*T^3 ...
                + 1.1407E-03*T^2 ...
                - 4.4890E-01*T ...
                + 1.0575E+03)/1000;                             % kJ/kg*K               // Air heat capacity
Cpxylene  =     (8.5210E4+3.2954E5*((1.4944E3/T)/sinh...
                (1.4944E3/T))^2+2.1150E5*((-6.7580E2/T)...
                /sinh(-6.7580E2/T))^2)/106.167E3;               % kJ/kg*K               // Xylene heat capacity

CpG       =     mair*Cpair+mxy*Cpxylene;                        % kJ/kg*K               // Gas heat capacity


% ----- EQUATIONS

%Subscript "ms" for the molten salt
%Subscript "g" for the gas

v0      =   2;                                                  % m/s                   // Molten salt mean velocity
Ab      =   pi*dti;                                             % m                     // Inside area (fluid's side)
Au      =   pi*dte;                                             % m                     // Outside area (molten salt's side)
Am      =   (Au-Ab)/log(Au/Ab);                                 % dimensionless         // Logaritmic mean between Au e Ab.
Gms     =   rhoms*us;                                           % kg/m²*s               // Massic velocity of molten salt
Gg      =   rhoG*v0;                                            % kg/m²*s               // Massic velocity of gas
Reg     =   (dp*Gg)/mig;                                        % dimensionless         // Reynolds number of gas
Rems    =   (dte*Gms)/mims;                                     % dimensionless         // Reynolds number of molten salt
Prms    =   (mims*cpms)/lambdams;                               % dimensionless         // Prandtl number of molten salt
epslon  =   0.38+0.073*(1+(((dt/dp)-2) ...
            ^2)/((dt/dp)^2));                                   % dimensionless         // Bed Porosity
f       =   6.8*((1-epslon)^1.2) ...
            *(Reg^-0.2)/(epslon^3);                             % dimensionless         // friction factor

        
% Defined dimensionless groups

A       =   0.62*Rems^0.5*Prms^(1/3);
B       =   (1+(0.4/Prms)^(2/3))^(1/4);
C       =   (1+(Rems/282000)^(5/8))^(4/5);
D       =   (3.50)*(Reg^0.7)*exp(-4.6*dp/dti);


alfau   =   (lambdams/dti)*(0.3+A/B*C);                         % W/m²*K                // Convective heat transfer coefficient from molten salt's side
alfai   =   (lambdag/dti)*D;                                    % W/m²*K                // Convective heat transfer coefficient from gas's side


Uinv    =   (1/alfai)+(d/lambdap)*(Ab/Am)+(1/alfau)*(Ab/Au);    %  1/U
U       =   (1/Uinv)/1000;                                      % kW/m²*K               // Global Heat Exchange Coefficient


% ---- END -----


