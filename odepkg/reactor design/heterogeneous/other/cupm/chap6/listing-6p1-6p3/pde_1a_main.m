%
% Wound healing
%
% Clear previous files
  clear all
  clc
%
% Global area
  global     nr     r0      r     u0...
              D      p     sc  ncall
%
% Model parameters
% r0=0.5; u0=1; D=2.0e-09; sc=      0; p=0; 
  r0=0.5; u0=1; D=2.0e-09; sc=8.0e-06; p=0; 
% r0=0.5; u0=1; D=2.0e-09; sc=8.0e-06; p=1; 
% r0=0.5; u0=1; D=2.0e-09; sc=8.0e-06; p=2;   
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
%
% Estimated velocity for p = 0
  for it=4:5
    fprintf('\n\n t = %6.2f',t(it));
    for i=1:nr
      fprintf('\n r = %5.3f   u(r,t) = %5.3f',r(i),u(it,i));
    end
  end 
  vs=1.0e+04*(0.175-0.120)/(1036800-777600);
  fprintf('\n\n velocity = %7.5f microns/s\n',vs);
  vday=(0.175-0.120)/(1036800-777600)*3600*24;
  fprintf('\n velocity = %6.4f cm/day\n',vday);
  rmon=(0.175-0.120)/(1036800-777600)*3600*24*30;
  fprintf('\n distance/month = %6.3f cm/mo\n',rmon);
  fprintf('\n wound half width = %6.3f cm\n',r0);
  fprintf('\n ncall = %5d\n',ncall);
  plot(r,u(2:nout,:),'-');
  axis([0 0.5 0 1])  
  xlabel('r (cm)'); ylabel('u(r,t)');
  title('u(r,t);t= 3, 6,..., 30 days') 
  
