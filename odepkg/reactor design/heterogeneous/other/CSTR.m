% CSTR.M
% Simple implementation of a transient Continuous Stirred Tank Reactor (CSTR) 
% carrying out a first-order reaction subjected to a first-order 
% catalyst deactivation, that is also dependent on the reactant concentration  

% Author: Luismar Porto
% (c) 2011 UFSC
% -------------------------------------------------------------------------
% The model is based on a molar balance for the key component A:
%
% dpsi/dtau = 1 - (1 + alpha*psi*a;
% 
% and a first order deactivation model
%
% da/dtau = -beta*psi*a;
%
% where:
% psi = dimensionless concentration
% a = catalyst activity
% alpha, beta => dimensionless parameters
%
% Initial conditions: psi = a = 1 for tau = 0

function CSTR           % Main function

clc;                    % Clear command window
clear all;              % Clear variables and functions from memory
close all; 
ic=[1; 1];              % Initial conditions
tspan = [0 2];          % Integration range

% Call the solver (ode45, ode15s, ...)

[t,x] = ode45 (@tank, tspan, ic);  

% Plot
hold on
figure (1)

% Plot concentration profile

subplot (1,2,1)
plot (t, x(:,1))
title ('Concentration profile')
xlabel ('Dimensionless time')
ylabel ('Dimensionless concentration')

% Plot conversion

subplot (1,2,2)
plot (t, x(:,2))
title ('Catalyst activity')
xlabel ('Dimensionless time')
ylabel ('Activity')

function xdot = tank (t, x) 

alpha = 2;
beta = 1;

xdot(1,:) = 1 - (1 + alpha*x(2))*x(1);

xdot(2,:) = -beta*x(2)*x(1);
