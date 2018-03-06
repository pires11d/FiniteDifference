% ========================================================================= 
% FEDERAL UNIVERSITY OF SANTA CATARINA - UFSC
% TECHNOLOGY CENTER - CTC
% CHEMICAL AND FOOD ENGINEERING DEPARTMENT - EQA
% INTEGRATED TECHNOLOGIES LABORATORY (InteLab), http://www.intelab.ufsc.br
% COURSE: EQA5409 CÁLCULO DE REATORES II
% CLASS: 0746 / 2010.1
% Professor LUISMAR MARQUES PORTO
% E-MAIL: luismar@intelab.ufsc.br
% ------------------------------------------------------------------------- 
% PROBLEM DEFINITION
%
% Non-isothermal spherical pellet with deactivation for 1st order reaction
% Pellet1ND.m
% 
% u(1) = concentration
% u(2) = temperature
%
% -------------------------------------------------------------------------
% REACTION: A --> Products
%                  
% REACTION PARAMETERS
% ra = -k*Ca        % Reaction rate
% DeltaH            % Heat of reaction
% 
% DESCRIPTION:
%
% The problem considers a 1st order reaction taking place in a spherical 
% catalyst particule, where diffusion resistance and deactivation is taken 
% into account. The pellet is also non-isothermal.
%
% Note: This program uses the PDEPE Partial Differential Equation solver
%
% PROCESS VARIABLES
%
% Cas   = catalyst pellet external surface concentration
% Ts    = catalyst pellet external temperature
% as    = catalyst pellet external surface activity
% Ca0   = catalyst pellet initial surface concentration
% T0    = catalyst pellet initial surface temperature
% a0    = catalyst pellet initial surface activity
% -------------------------------------------------------------------------
% VARIABLES
%
% Ca    = molar concentration of component A (dependent)
% T     = temperature (dependent)
% a     = catalyst activity (dependent)
% t     = time (independent)
% x     = radius (independent)
% -------------------------------------------------------------------------
% Catalyst parameters
%
% k     = 1st order volumetric rate constant
% k0    = kinetic rate pre-exponential factor (Arrhenius constant)
% E     = main reaction (conversion of A) activation energy
% -------------------------------------------------------------------------
% Catalyst porous medium parameters
% De    = effective diffusivity in the catalyst pores
% ke    = effective thermal conductivity in the catalyst pores
% alphae = effective thermal difusivity in the catalyst pores
% -------------------------------------------------------------------------
% Fluid parameters
%
% rho   = molar specific mass (molar density)
% cp    = molar heat capacity
% -------------------------------------------------------------------------
% Catalyst pellet parameters
%
% Rc    = catalyst radius
% De    = catalyst effective diffusivity inside the pores
% -------------------------------------------------------------------------
% Deactivation parameters
% 
% kd    = 1st order volumetric deactivation rate constant
% kd0   = kinetic deactivation rate pre-exponential factor
% Ed    = activation energy of the deactivation process
% -------------------------------------------------------------------------
% DIMENSIONLESS GROUPS
%
% phi   = (Rc^2*k0*exp(E/R*Ts))^(1/2)   % Thiele modulus for the main rxn
% phia  = (Rc^2*kd0*exp(Ed/R*Ts))^(1/2) % Thiele modulus for deactivation
% beta  = -DeltaH*De*Cas/(k*Ts)         % Temperature elevation number
% gamma = E/(R*Ts)                      % Arrhenius number
% Le    = (ke/(rho*cp))/De = alphae/De  % Modified Lewis number
% -------------------------------------------------------------------------
% Dimensionless quantities and definitions
%
% psi   = Ca/Cas;               % Dimensionless concentration
% theta = T/Ts;                 % Dimensionless temperature
% a     = (-ra)/(-ra0)          % Catalyst activity
% qsi   = r/R                   % Pellet dimensionless (relative) radius
% tau   = t*De/R^2              % Dimensionless time
% -------------------------------------------------------------------------
% Mole balance for component A:
% 
% Dpsi/Dtau = (1/qsi^2)*D/Dqsi*(qsi^2*Dpsi/Dqsi) - ...
% theta^2*psi*exp[(1-1/theta)]*exp(-phia^2*tau)
% 
% Energy (heat) balance
% 
% (1/Le)*Dtheta/Dtau = (1/qsi^2)*Dtheta/Dqsi*(qsi^2*Dtheta/Dqsi) + ...
% beta*theta^2*psi*exp((gamma(1-1/theta))*exp(-phia^2*tau)
% 
% Initial conditions:
% psi(0,qsi) = 1, theta(0,qsi) = 1; a(0,qsi) = 1
%
% Boundary conditions:
% psi(tau,1) = 1, theta(tau,1) = 1; |a(tau,1) = 1|
% Dpsi/Dqsi(tau,0) = 0, Dtheta/Dqsi(tau,0) = 0, |Da/Dqsi(tau,0) = 0|
%
% Main function
%
function Pellet1NDmod % 1st order, Non-isothermal, with Deactivation
close;
m = 2; % spherical coordinates
x = linspace(0,1,100); % We are taking 100 time points, from 0 to 1
t = linspace(0,1,100);  % We are taking  20 space points (radial), from 0 to 1

% Call the PDE solver, passing the equations (coefficients, initial and
% boundary conditions
%
sol = pdepe(m,@Pellet1NDpde,@Pellet1NDic,@Pellet1NDbc,x,t);
% Extract the u components
u1 = sol(:,:,1);    % concentration
u2 = sol(:,:,2);    % temperature
u3 = sol(:,:,3);    % activity

% Concentration

surf(x,t,u1)    
title('Dimensionless concentration')
xlabel('Distance r')
ylabel('Time t')
% A solution profile can also be illuminating

figure (2)
plot(x,u1(end,:))
title('Solution at the final t')
xlabel('Distance r')
ylabel('c')

% Temperature
surf(x,t,u2)    
title('Dimensionless temperature')
xlabel('Distance r')
ylabel('Time t')
% A solution profile can also be illuminating

figure (3)
plot(x,u2(end,:))
title('Solution at the final t')
xlabel('Distance r')
ylabel('T')
% -------------------------------------------------------------------------
% 
function [c,f,s] = Pellet1NDpde(x,t,u,DuDx)

% Parameters for the synthesis of vinyl chloride
phi = 0.27; beta = 0.25; gamma = 6.5; Le = 0.1; phia = 0.1*phi;

% Parameters for the synthesis of vinyl chloride
%phi = 5.0; beta = 0.64; gamma = 22.0; Le = 0.01; phia = 1*phi;

% Cannonical coefficients

c = [1; 1/Le; 1];
f = [1; 1; 0].*DuDx;
s = [-phi^2.*u(1).*exp(1-1/u(2)).*u(3); ...
    beta*phi^2.*u(1)*exp(gamma*(1-1/u(2))).*u(3); ...
    -phia^2.*u(3)];
% -------------------------------------------------------------------------
function u0 = Pellet1NDic(x)
u0 = [0; 1; 1]; 
% -------------------------------------------------------------------------
function [pl,ql,pr,qr] = Pellet1NDbc(xl,ul,xr,ur,t)
pl = [0; 0; ul(3)-1];
ql = [1; 1; 0];
pr = [ur(1) - 1; ur(2) - 1; ur(3)-1]; 
qr = [0; 0; 0];  