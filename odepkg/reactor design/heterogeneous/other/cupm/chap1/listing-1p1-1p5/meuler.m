  function [u,t]=meuler(u0,t0,h,nsteps)
%
% nsteps modified Euler steps
  for i=1:nsteps
%
%   Modified Euler integration
%     ut=pde_1(u0,t0); 
%     u1=u0+ut*h;
%     t=t0+h;
%     u1t=pde_1(u1,t);
%     u=u0+(ut+u1t)*h/2;
%
%   Runge Kutta format
      k1=pde_1(u0,t0); 
      u1=u0+k1*h;
      t=t0+h;
      k2=pde_1(u1,t);
      u=u0+(k1+k2)*h/2;
%
% Next modified Euler step
  u0=u; t0=t;
  end  
