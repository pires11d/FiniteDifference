%
% Clear previous files
  clear all
  clc
%
% Parameters shared with the ODE routine
  global zl zu z dz dz2 D kf cbsat kr n cbulk ndss ncall
%
% Parameter numerical values
  D=1.0e-10; kf=1.0e+05; cbulk=4.48e-05; 
  h=5.0e-05; c0=0; cb0=0;
%
% Variation in interface binding saturation
  cbsat=1.66e-08;
  cbsat=1.66e-09;
%
% Variation in interface unbinding rate
  kr=1.0e-01;
  kr=1.0e+01;
%
% Spatial grid
  zl=0; zu=5.0e-05; n=21; dz=(zu-zl)/(n-1); dz2=dz^2;
%
% Initial condition
  for i=1:n
    u0(i)=c0; 
  end 
  u0(n+1)=cb0;
%
% Independent variable for ODE integration
  t0=0.0;
  tf=100;
  tout=(t0:2:tf); 
  nout=51;
  ncall=0;
%
% ODE itegration
%
% Variation in error tolerances
  reltol=1.0e-06; abstol=1.0e-06;
  reltol=1.0e-07; abstol=1.0e-07;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
  mf=1;
  if(mf==1) % explicit FDs 
    [t,u]=ode15s(@pde_1,tout,u0,options); end 
  if(mf==2) ndss=4; % ndss = 2, 4, 6, 8 or 10 required
    [t,u]=ode15s(@pde_2,tout,u0,options); end
  if(mf==3) ndss=44; % ndss = 42, 44, 46, 48 or 50 required 
    [t,u]=ode15s(@pde_3,tout,u0,options); end
%       
% Store numerical solutions
  for it=1:nout
    c_plot(it)=u(it,1); 
    cb_plot(it)=u(it,n+1);
    theta_plot(it)=cb_plot(it)/cbsat;
    rate_plot(it)=kf*c_plot(it)*(cbsat-cb_plot(it))-kr*cb_plot(it);
  end
%
% Display selected output
  fprintf('\n mf = %2d   abstol = %8.1e   reltol = %8.1e\n',...
          mf,abstol,reltol);
  fprintf('\n     t      c(0,t)       cb(t)       theta        rate\n');
  for it=1:nout
    fprintf('%6.0f%12.3e%12.3e%12.3e%12.3e\n',...
            t(it),c_plot(it),cb_plot(it),theta_plot(it),rate_plot(it));
  end
  fprintf('\n ncall = %4d\n',ncall);
%
% Plot numerical solution
  figure(1);
  subplot(2,2,1)
  plot(t,c_plot); axis tight
  title('c(0,t) vs t'); xlabel('t'); ylabel('c(0,t)')
  subplot(2,2,2)
  plot(t,cb_plot); axis tight
  title('cb(t) vs t'); xlabel('t'); ylabel('cb(t)')
  subplot(2,2,3)
  plot(t,theta_plot); axis tight
  title('theta(t) vs t'); xlabel('t'); ylabel('theta(t)')
  subplot(2,2,4)
  plot(t,rate_plot); axis tight
  title('rate(t) vs t'); xlabel('t'); ylabel('rate(t)')
%  
% Store numerical solution for 3D plot
  for it=1:nout
  for i=1:n    
    c_3D(it,i)=u(it,i); 
  end
  end
  z=[zl:dz:zu];
  figure(2)  
  surf(z,t,c_3D)
  xlabel('z (m)'); ylabel('t (sec)'); zlabel('c(z,t) (moles/m^3)');
  title('c(z,t) (moles/m^3), z=0,2.5\times10^{-6},..., 5\times10^{-5} (m), t=0,2,...,100 (sec)')
% print -deps pde.eps; print -dps pde.ps
