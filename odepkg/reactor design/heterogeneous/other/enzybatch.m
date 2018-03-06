%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% UNIVERSIDADE FEDERAL DE SANTA CATARINA - UFSC
% CENTRO TECNOLÓGICO - CTC
% DEPARTAMENTO DE ENGENHARIA QUÍMICA E ENGENHARIA DE ALIMENTOS - EQA
% DISCIPLINA: EQA5409 - CÁLCULO DE REATORES II
% TURMA: 0746 / 2009.2
% PROFESSORES:                  LUISMAR MARQUES PORTO & CÍNTIA SOARES
% ALUNA DE ESTÁGIO DOCÊNCIA:    FERNANDA VIEIRA BERTI
% E-MAIL: luismar@intelab.ufsc.br
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PROBLEM DEFINITION
%
% ENZYBATCH.m
%
% BATCH REACTOR WITH ENZYMATIC REACTION
%
% REACTION: A --> Products
%                  
% x(1)= psi     (dimensionless concentration)
% x(2)= a       (dimensionless activity)
%
% 
% xdot(1)= dpsi/dtheta
% xdot(2)= da/dtheta
% 
% DESCRIPTION:
% The problem considers a batch reactor of fixed (liquid) volume, V, where
% an enzymatic hydrolysis reaction takes place. The enzyme activity, for 
% some undefined reason, decays according to a deactivation rate that is 
% dependent on the activity itself and the substrate concentration.
% The main reaction follows a Langmuir-Hinshelwood-Hougen-Watson (LHHW)
% mechanism, with a relatively weak interaction (adsorption) constant.
% Below are the two major ODEs that govern the batch reaction system, for a
% isothermal reactor.
%
% Mole balance:
% dpsi/dtheta = - alpha*psi*a/(1 + beta*psi)
% 
% Activity equation:
% da/dtheta = - gamma*psi*a
% 
% For theta = 0, psi = 1, a = 1
%
% The model assumes dimensionless concentration and enzyme activity.
% alpha = k*tf*Ca0/V
% beta  = Ka*Ca0
% gamma = tf*Ca0
% 
% where:
%           k   = kinetic constant, L/h
%           tf  = batch time, h
%           Ca0 = feed concentration, mol/L
%           V   = batch volume, L
%           Ka  = adsorption constant, L/mol
%
% ----- Part 1: Main function, cleaning up, initial conditions, call to 
% the ODE solver

function enzybatch      % Main function

clc;                    % Clear command window
clear;                  % Clear variables and functions from memory
close all;              % Closes all the open figure windows

% Initial conditions and integration range

ic=[1;1];               % Initial conditions
thetaspan = [0 1];      % Integration range

% Call the solver (ode45, ode15s, ...)

[theta,x] = ode45(@enzyme, thetaspan, ic);  

% Plot

hold on
figure (1)

% Plot concentration profile
plot(theta, x(:,1),'b:o', theta, x(:,2), 'r:o');
title ('Concentration and Activity Profile')
xlabel ('Batch time, Dimensionless')
ylabel ('Concentration and Activity, Dimensionless')  
legend('Concentration', 'Activity','location', 'southwest')
grid

% ----- Part 2: Problem definition

function xdot = enzyme(theta, x) 

% Dimensionless parameters

alpha =1;
beta = 0.01;
gamma =0.1;

% ----- Part 3: Set of dimensionless ordinary diffential equations (ODEs)

xdot(1, :)= -alpha*x(1)*x(2)/(1+beta*x(1));    % Mass (or molar) balance
xdot(2, :)= -gamma*x(1)*x(2);                  % Activity

% ---- END -----