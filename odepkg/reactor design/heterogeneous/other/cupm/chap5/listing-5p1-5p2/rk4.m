  function [u,t]=rk4(u0,t0,h,nsteps)
%
% nsteps modified Euler steps
  for i=1:nsteps
%
%   Runge Kutta integration
      k1=pde_1(u0,t0);
      u1=u0+k1*h/2;
      t=t0+h/2;
      k2=pde_1(u1,t);
      u1=u0+k2*h/2;
      k3=pde_1(u1,t);
      u1=u0+k3*h;
      t=t0+h;
      k4=pde_1(u1,t);      
      u=u0+(1/6)*(k1+2*k2+2*k3+k4)*h;
%
% Next Runge Kutta step
  u0=u; t0=t;
  end  
