%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% UNIVERSIDADE FEDERAL DE SANTA CATARINA - UFSC
% CENTRO TECNOLÓGICO - CTC
% DEPARTAMENTO DE ENGENHARIA QUÍMICA E ENGENHARIA DE ALIMENTOS - EQA
% DISCIPLINA: EQA5409 - CÁLCULO DE REATORES II
% TURMA: 0746 / 2010.1
% PROFESSOR LUISMAR MARQUES PORTO
% E-MAIL: luismar@intelab.ufsc.br
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PROBLEM DEFINITION
%
% PelletO1st.m
%
% SPHERICAL PELLET WITH DIFFUSION RESISTANCE AND A 1ST ORDER REACTION
%
% REACTION: A --> Products
%                  
% u = concentration
%
% DESCRIPTION:
% The problem considers a 1st order reaction taking place in a spherical 
% catalyst particule, where diffusion resistance is taken into account
%
% Mole balance:
% 
% DCa/Dt + div(-D*grad(Ca)) -k*Ca = 0
% 
% Initial condition:
% Ca = 0 at t = 0
%
% Boundary conditions:
% DCa/Dt = 0 at r = 0
% Ca = Ca0 at r = R
%
% Main function
%
function pelletO1st
close;
m = 2; % spherical coordinates
x = linspace(0,1,100); % We are taking 100 time points, from 0 to 1
t = linspace(0,1,20);  % We are taking  20 space points (radial), from 0 to 1

% Call the PDE solver, passing the equations (coefficients, initial and
% boundary conditions
%
sol = pdepe(m,@pelletO1pde,@pelletO1ic,@pelletO1bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);

% A surface plot is often a good way to study a solution.
surf(x,t,u)    
title('Numerical solution computed with 20 mesh points.')
xlabel('Distance r')
ylabel('Time t')
% A solution profile can also be illuminating.

figure (2)
plot(x,u(end,:))
title('Solution at t = 1')
xlabel('Distance r')
ylabel('u(x,2)')


% --------------------------------------------------------------
function [c,f,s] = pelletO1pde(x,t,u,DuDx)
D = 1;
k = 1;
c = 1;
f = D*DuDx;
s = -k*u;
% --------------------------------------------------------------
function u0 = pelletO1ic(x)
u0 = 0;
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pelletO1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = -1;
qr = 1; 