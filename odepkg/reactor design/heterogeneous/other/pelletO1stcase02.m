%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% UNIVERSIDADE FEDERAL DE SANTA CATARINA - UFSC
% CENTRO TECNOL?GICO - CTC
% DEPARTAMENTO DE ENGENHARIA QU?MICA E ENGENHARIA DE ALIMENTOS - EQA
% DISCIPLINA: EQA5409 - C?LCULO DE REATORES II
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
% DCa/Dr = 0 at r = 0
% Ca = Ca0 at r = R
%
% Main function
%
function pelletO1stcase02
close;
m = 2; % spherical coordinates
x = linspace(0,1,100); % We are taking 100 time points, from 0 to 1
t = linspace(0,1,20);  % We are taking  20 space points (radial), from 0 to 1

% Call the PDE solver, passing the equations (coefficients, initial and
% boundary conditions
%
sol = pdepe(m,@pelletO1pde,@pelletO1ic,@pelletO1bc,x,t);
% Extract the first solution component as u.
u1 = sol(:,:,1);
u2 = sol(:,:,2);
u3 = sol(:,:,3);


% Concentration.

surf(x,t,u1)    
title('Numerical solution computed with 20 mesh points.')
xlabel('Distance r')
ylabel('Time t')
% A solution profile can also be illuminating.

figure (2)
plot(x,u1(end,:))
title('Solution at t = 1')
xlabel('Distance r')
ylabel('u1(x,1)')

% Temperature.
surf(x,t,u2)    
title('Numerical solution computed with 20 mesh points.')
xlabel('Distance r')
ylabel('Time t')
% A solution profile can also be illuminating.

figure (4)
plot(x,u2(end,:))
title('Solution at t = 1')
xlabel('Distance r')
ylabel('u2(x,1)')

% Activity.
surf(x,t,u3)    
title('Numerical solution computed with 20 mesh points.')
xlabel('Distance r')
ylabel('Time t')
% A solution profile can also be illuminating.

figure (6)
plot(x,u3(end,:))
title('Solution at t = 1')
xlabel('Distance r')
ylabel('u3(x,1)')



% --------------------------------------------------------------
function [c,f,s] = pelletO1pde(x,t,u,DuDx)
Cas = 1; % not used
Ts = 1;  % not used
a = 1;   % not used
D = 1;
lambda = 1;
fakeD = 0; % fake diffusivity for activity
k = 1;
kd = 1;
rho = 1;
cp = 1;
DeltaH = -1;

c = [1; rho*cp; 1];
f = [D; lambda; fakeD].*DuDx;
s = [-k*u(1)*u(3); -DeltaH*k*u(1)*u(3); -kd*u(3)]; 
% --------------------------------------------------------------
function u0 = pelletO1ic(x)
u0 = [0; 1; 1]; 
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pelletO1bc(xl,ul,xr,ur,t)
pl = [0; 0; 0];
ql = [1; 1; 1];
pr = [ur(1) - 1; ur(2) - 1; ur(3)]; 
qr = [0; 0; 0];  