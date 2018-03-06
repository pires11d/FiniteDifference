%
% Wound healing
%
% Clear previous files
  clear all
  clc
  global     nr     r0      r     u0...
              D      p     sc  ncall
%
% Model parameters
  r0=0.5; u0=1; D=2.0e-09; sc=      0; p=0; 
% r0=0.5; u0=1; D=2.0e-09; sc=8.0e-06; p=0; 
% r0=0.5; u0=1; D=2.0e-09; sc=8.0e-06; p=1; 
% r0=0.5; u0=1; D=2.0e-09; sc=8.0e-06; p=5;   
%
% Spatial grid
  nr=101;
  dr=r0/(nr-1); 
  for j=1:nr
    r(j)=(j-1)*dr;
  end
%
% Independent variable for ODE integration
  tf=60*60*24*30;
  tout=[0:tf/10:tf]'; 
  nout=11;
  ncall=0;
%
% Initial condition
  for i=1:nr u0w(i)=0; end
%
% ODE integration 
  reltol=1.0e-04; abstol=1.0e-04;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
  [t,u]=ode15s(@pde_1,tout,u0w,options);   
%
% Display a heading and t at output
  fprintf('\n  nr = %2d\n',nr);
  for it=1:nout
    u(it,1)=u0;
    fprintf('\n t (s) = %6.2e   t (days) = %6.2f',...
            t(it),t(it)/(60*60*24));
  end
  fprintf('\n\n ncall = %5d\n',ncall);
%
% Plot the solution as u(r,t) vs r with t as a parameter
  plot(r,u(2:nout,:),'-');
  axis([0 0.5 0 1])
  xlabel('r (cm)'); ylabel('u(r,t)');
  title('u(r,t);t = 3, 6,..., 30 days') 