  clc
  clear all
%
%  Linear advection equation
%
%  The linear advection equation
%
%    ut + v*uz = 0                                          (1)
%
%  is integrated by the method of lines (MOL) subject to 
%  the IC
%
%    u(z,t=0) = f(z)                                        (2)
%  
%  BC
%
%    u(z=0,t) = g(t)                                        (3)
%
%  We consider in particular f(z) = 0, g(t) = 1 corresponding 
%  to the Heaviside unit step function, h(t); the solution to 
%  eqs. (1) to (3) is
%
%    u(z,t)=h(t - z/v)                                      (4)
%
% which is used to evaluate the numerical (MOL) solution.
%
%  The numerical algorithms are:
%
%    z (spatial, boundary value) integration: Two point upwind
%      (2pu)
%
%    t (temporal, initial value) integration: Explicit Euler
% 
  global z dz zL v ifd n ncase ncall
%
% Select an approximation for the convective derivative uz
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
% Grid (in z)
  zL=1;  n=51; dz=0.02;
% zL=1; n=101; dz=0.01;  
  z =[0:dz:zL];
%
% Step through cases
%
%   ncase = 1: v > 0
%
%   ncase = 2: v < 0
%
  for ncase=1:2
    if ncase==1 v= 1; end
    if ncase==2 v=-1; end
%
% Parameters for Euler integration
  nsteps=250;
  h=0.001;
%
% Initial condition
  if(ncase==1)
    for i=2:n
      u(i)=0;
    end
    u(1)=1;
  end 
  if(ncase==2)
    for i=1:n-1
      u(i)=0;
    end
    u(n)=1;
  end    
  t=0;
%
% Write ncase, h, v
  fprintf('\n\n ncase = %5d    h = %10.3e   v = %4.2f\n\n',ncase,h,v);
%
% nout output points
  nout=3;
  ncall=0;
  for iout=1:nout
%
%   Euler integration
    u0=u; t0=t;
    [u,t]=euler(u0,t0,h,nsteps);
%
%   Store solution for plotting at t = 0.25, 0.5, 0.75
    if(ncase==1)
      for i=1:n
        if(z(i)>= v*t) ua(i)=0;   end
        if(z(i)<  v*t) ua(i)=1;   end
         uplot(iout,i)=u(i);
        uaplot(iout,i)=ua(i);      
      end  
    end 
    if(ncase==2)
      for i=1:n
        if(z(i)>=z(n)-abs(v)*t) ua(i)=0;   end
        if(z(i)< z(n)-abs(v)*t) ua(i)=1;   end
         uplot(iout,i)=u(i);
        uaplot(iout,i)=ua(i);      
      end  
    end 
%
% Next output
  end
%
% Plots for ncase = 1, 2
  figure(ncase);
  plot(z,uplot,'-o');
  if(ifd==1)axis([0 1 0 1]); end
  if(ifd==2)axis([0 1 0 1.5]); end
  if(ifd==3)axis([0 1 -0.5 1.5]); end
  if(ifd==4)axis([0 1 0 1]); end   
  if(ifd==5)axis([0 1 0 1]); end  
  if(ifd==6)axis([0 1 0 1]); end      
  ylabel('u(z,t),ua(z,t)');xlabel('z');
  if(ncase==1)
    title('ncase = 1, v > 0; t = 0.25, 0.5, 0.75; num - o; anal - line');
  end
  if(ncase==2)
    title('ncase = 2, v < 0; t = 0.25, 0.5, 0.75; num - o; anal - line');
  end
  hold on
  plot(z,uaplot,'-');
%
% Next case
  end
