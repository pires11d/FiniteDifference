  function [u,t]=euler(u0,t0,h,nsteps)
%
% nsteps Euler steps
  for i=1:nsteps
%
%   Euler integration
%     ut=pde_1(u0,t0); 
%     u=u0+ut*h;
%     t=t0+h;      
%
%   Runge Kutta format - 1
%     k1=pde_1(u0,t0); 
%     u=u0+k1*h;
%     t=t0+h;      
%
%   Runge Kutta format - 2
      k1=pde_1(u0,t0)*h; 
      u=u0+k1;
      t=t0+h;      
%
% Next Euler step
  u0=u; t0=t;
  end  