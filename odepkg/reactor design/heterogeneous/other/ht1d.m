% HT1D.M
% Heterogeneous 1D Reactor Model
function ht1d

close all;
clc;



% Integrate concentration profile to obtain effectiveness factor

m = 0; % axial flow with dispersion
x_initial   =  1;
x_final     = 10;
t_initial   =  0;
t_final     =  2;


mesh_x = 100;   % Number of mesh points for x
mesh_t = 3;    % Number of mesh point for t
delta_t = (t_final - t_initial)/mesh_t; % Delta t

x = linspace(x_initial,x_final,mesh_x);  % We are taking 100 time points, from 1 to 20
%t_initial = time;

%t = linspace(t_initial,t_final,mesh_t);  % We are taking  50 space points (radial), from 0 to 2
% Call particle balances
% Do we need a call for first particle?? (Do not forget to define t1 and t2
% particle = f_particle (0, 1); % particle contains c, T and a profiles

% Call the PDE solver, passing the equations (coefficients, initial and
% boundary conditions
%

time = 0;
for i = t_initial : t_final : mesh_t

    t1 = t_initial;
    t2 = t1 + delta_t;
    t = linspace(t1,t2,mesh_t);
    
    reactor = pdepe(m,@ht1dpde,@ht1dic,@ht1dbc,x,t);

    % Extract the first solution component as u.
    u = reactor(:,:,1);
    particle = f_particle (t1, t2);
    time = time + delta_t;
end
% A surface plot is often a good way to study a solution.
surf(x,t,u)    
title('Numerical solution computed with 100 mesh points.')
xlabel('Distance r')
ylabel('Time t')
% A solution profile can also be illuminating.

figure (2)
plot(x,u(end,:))
title('Solution at t = 1')
xlabel('Distance r')
ylabel('psi')

% --------------------------------------------------------------
function [c,f,s] = ht1dpde(x,t,u,DuDx)

% Parameters
Peam = f_Peam;
Peah = f_Peah;
eta = f_eta;
rib = f_rib;
DeltaHj = f_DeltaHj;
delta = f_delta;
etaj = f_etaj;

% ------------  coefficients ----

c = [ 1;                    1;          0];
f = [(1/Peam);              (1/Peah);   0].*DuDx;
s = [ eta*rib*u(3) - DuDx;   DeltaHj*etaj*rib*u(3) - DuDx; -DuDx + delta ];
% --------------------------------------------------------------
function u0 = ht1dic(x)
u0 = [1;1;1];
% --------------------------------------------------------------
function [pl,ql,pr,qr] = ht1dbc(xl,ul,xr,ur,t)

pl = [1;1;1];
ql = [1;1;1];
pr = [1;1;1];
qr = [1;1;1]; 

% ==============================================================

function Peam = f_Peam
...
Peam = 1;

function Peah = f_Peah
...
Peah = 1;
function eta = f_eta
...
eta = 1;
function etaj = f_etaj
...
etaj = 1;
function rib = f_rib
...
    rib = 1;
function DeltaHj = f_DeltaHj
...
    DeltaHj = 1;
function delta = f_delta
...
    delta = 1;
function [c, T, a] = f_particle (time1, time2)

m = 2; % spherical particle with radial symmetry
mesh_rp = 100;   % Number of mesh points for rp
mesh_tp = 50;    % Number of mesh points for tp

rp = linspace(1,10,mesh_rp); %

time1 = time1; % Need correction
time2 = time2; % Need correction

tp = linspace(time1,time2,mesh_tp);  % 

% Call the PDE solver (pdepe), passing the equations (coefficients, initial and
% boundary conditions)
%
pellet = pdepe(m,@pelletpde,@pelletic,@pelletbc,rp,tp);
% Extract the first solution component as c, second as T and third as a.
c = pellet(:,:,1);
T = pellet(:,:,2);
a = pellet(:,:,3);

% A surface plot is often a good way to study a solution.
%surf(x,t,c)    
%title('Numerical solution computed with 100 mesh points.')
%xlabel('Distance r')
%ylabel('Time t')
% A solution profile can also be illuminating.

%figure (2)
%plot(x,c(end,:))
%title('Solution at t = 1')
%xlabel('Distance r')
%ylabel('psi')
function [cp,fp,sp] = pelletpde(rp,tp,c,DcDx)

% ------------  coefficients ----
cp = [1 ; 1; 1];
fp = [1 ; 1; 1].*DcDx;
sp = [1 ; 1; 1];
% --------------------------------------------------------------
function c0 = pelletic(r)
c0 = [1;1;1];
% --------------------------------------------------------------
function [plp,qlp,prp,qrp] = pelletbc(rl,cl,rr,cr,tp)

plp = [1;1;1];
qlp = [1;1;1];
prp = [1;1;1];
qrp = [1;1;1]; 
