  clc
  clear all
%
%   1D dialyzer model
%
%   The ODE/PDE system is 
%
%   u1_t = -v1*u1_z + kM1*(u2 - u1)                        (PDE)(1a)
%
%   u2_t = -v2*u2_z + kM2*(u1 - u2)                        (PDE)(1b)
%
%   V1L*u1_t(z=0,t) = q1*(u1L - u1(z=0,t))                 (ODE)(1c)
%
%   V1R*u1R_t(t) = q1*(u1(z=zL,t) - u1R(t))                (ODE)(1d)
%
%   For the countercurrent case (v1 > 0, v2 < 0).
%
%   The primary outputs are u1R(t), u2(z=0,t) as a function t.
%
%   The method of lines (MOL) solution for eqs. (1) is coded below.  
%   Specifically, the spatial derivatives in the fluid balances, u1_z,
%   u2_z, are replaced by one of six approximations as selected by the 
%   variable ifd.
%
%   The following cases are programmed:
%
%   ncase = 1: kM = 0 (no mass transfer)
%
%   ncase = 2: kM = 0.001 (blood to dialysate mass transfer)
%
  global zL q1 q2 v1 v2 V1L V1R kM1 kM2 u1L u2zL u1 u2 ifd n ncase ncall
%
% Step through cases
  for ncase=1:2
%
%   Model parameters
    D=5; A=pi*D^2/4; AM=A; q1=0.25; q2=0.25; eps=0.5; 
    v1=q1/(eps*A); v2=-q2/((1-eps)*A); u1L=1; 
    V1L=A; V1R=A; u10=0; u20=0; zL=50; n=21;
%
%   Set parameters for each case
    if(ncase==1) kM=0;      u2zL=1; end
    if(ncase==2) kM=0.001;  u2zL=0; end
    kM1=AM*kM/(A*eps); 
    kM2=AM*kM/(A*(1-eps));
%
%   Display parameter summary
    fprintf(...
    '\n\n D = %3.1f  A = %4.1f  q1 = %4.2f  q2 = %4.2f  eps = %4.2f\n',...
    D,A,q1,q2,eps);
    fprintf(...
    '\n v1 = %5.3f  v2 = %5.3f  u1L = %5.3f  u2zL = %4.2f  kM = %5.3f\n',...
    v1,v2,u1L,u2zL,kM);
    fprintf(...
    '\n u10 = %4.2f  u20 = %4.2f  zL = %4.1f  n = %3d\n',...
    u10,u20,zL,n);    
%
% Select an approximation for the convective derivatives u1z, u2z
%
%   ifd = 1: Two point upwind approximation
%
%   ifd = 2: Centered approximation
%
%   ifd = 3: Five point, biased upwind approximation
%
%   ifd = 4: van Leer flux limiter
%
%   ifd = 5: Superbee flux limiter
%
%   ifd = 6: Smart flux limiter
%
  ifd=1;
%
% Level of output
%
%   Detailed output - ip = 1
%
%   Brief (IC) output - ip = 2
%
  ip=1;
%
% Parameters for fourth order Runge Kutta integration
  nsteps=10;
  h=14.4;
%
% Initial condition
  for i=1:n
      u(i)=u10;
    u(i+n)=u20;
  end
  u(2*n+1)=u10;
  t=0;
%
% Display ifd, ncase, h, CFL
  fprintf('\n ifd = %2d   ncase = %2d   h = %10.3e  CFL = %4.2f\n\n',...
          ifd,ncase,h,v1*h/(zL/(n-1)));
%
% Display heading
  fprintf('    t    u1R(t)   u2(0,t)\n');
%
% Display numerical solution at t = 0
  fprintf('%5.2f%10.4f%10.4f\n',t/3600,u(2*n+1),u(n+1));
%
% Store solution for plotting
   u1plot(1)=u(2*n+1);
   u2plot(1)=u(n+1);
   tplot(1)=t;
%
% nout output points
  nout=51;
  ncall=0;
  for iout=2:nout
%
%   Fourth order Runge Kutta integration
    u0=u; t0=t;
    [u,t]=rk4(u0,t0,h,nsteps);      
%
%   Numerical solutions
    if(ip==1)
      fprintf('%5.2f%10.4f%10.4f\n',t/3600,u(2*n+1),u(n+1));
    end  
%
%   Store solution for plotting
     u1plot(iout)=u(2*n+1);
     u2plot(iout)=u(n+1);
     tplot(iout)=t/3600;    
%
% Next output
  end
%
% Calls to ODE routine
  fprintf('\n ncall = %4d\n',ncall);
%
% Plots for u1R(t), u2(z=0,t)
  figure(ncase);
  plot(tplot,u1plot,'-o');
  axis([0 2 -0.1 1.1]);
  ylabel('u1R(t),u2(0,t)');xlabel('t');
  if(ncase==1) title('ncase = 1; u1 - o; u2 - x'); end
  if(ncase==2) title('ncase = 2; u1 - o; u2 - x'); end
  hold on
  plot(tplot,u2plot,'-x');
%
% Next case
  end
