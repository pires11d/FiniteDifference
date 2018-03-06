  clc
  clear all
% 
%   One dimensional convection diffusion reaction (CDR) equation
%
%   The PDE is
%
%      ut = D*uzz - v*uz - kr*u                             (1)
%
%   with the initial condition
%
%     u(z,t=0) = f(z)                                       (2)
%
%   and boundary conditions
%
%     v*u(z=0,t) = v*ue - D*uz(z=0,t);                      (3)
%
%     uz(z=zl,t) = 0                                        (4)
%
%   where
%
%      u     dependent variable, e.g., concentration
%
%      t     time
%
%      z     position
%
%      D     diffusivity (dispersion coefficient)
%
%      v     velocity
%
%      kr    reaction rate constant
%    
%      f(z)  initial condition function
%    
%      ue    entering u
%
%      zl    length
%
%   The following code is run for five special cases:
%
%      Case 1: uz by two point upwind (2pu), uzz by three point 
%              centerd finite differences (3pc) and direct second 
%              order differentiation
%
%      Case 2: uz and uzz by three point centered finite differences 
%              (3pc) and direct second order differentiation
%
%      Case 3: uz by van Leer flux limiter (vanl), uzz by three 
%              point centered finite differences (3pc) and direct 
%              second order differentiation
%
%      Case 4: uz by superbee flux limiter (super), uzz by five 
%              point centered finite differences (5pc) and direct 
%              second order differentiation
%
%      Case 5: uz by smart flux limiter (smart), uzz by seven 
%              point centered finite differences (7pc) and direct 
%              second order differentiation
%
  global a b D v ue kr zl ncase n ncall
%
% Model parameters (defined in the comments above)
  kr=0; v=1; ue=1; zl=1; n=51; D=0.0001;
%
% Step through cases
  for ncase=1:5
%
% Level of output
%
%   Detailed output - ip = 1
%
%   Brief (IC) output - ip = 2
%
  ip=2;
%
% Parameters for fourth order Runge Kutta integration
  nsteps=200;
  h=0.001;
%
% Initial condition
  for i=1:n
    z(i)=(i-1)/(n-1)*zl;  
    u(i)=0;
  end
  u(1)=ue;
  t=0;
%
% Write ncase, h
  fprintf('\n\n ncase = %2d   h = %10.3e\n\n',ncase,h);
%
% nout output points
  nout=6;
  ncall=0;
  for iout=1:nout
%
%   Fourth order Runge Kutta integration
    u0=u; t0=t;
    [u,t]=rk4(u0,t0,h,nsteps);      
%
%   Numerical solution
    fprintf('\n t = %5.2f\n',t);
    for i=1:n
      uplot(iout,i)=u(i);
      if(ip==1)
        fprintf('%7.2f%12.3f\n',z(i),u(i));
      end  
    end 
%
% Next output
  end
%
% Plots for ncase = 1 to 5
  figure(ncase);
  plot(z,uplot,'-o');
  axis([0 1 0 1.2]);
  ylabel('u(z,t)'); xlabel('z');
  if(ncase==1)title('ncase = 1, u(z,t), t = 0.2, 0.4, ..., 1.2,');end
  if(ncase==2)title('ncase = 2, u(z,t), t = 0.2, 0.4, ..., 1.2,');
                                                axis([0 1 0 1.4]);end
  if(ncase==3)title('ncase = 3, u(z,t), t = 0.2, 0.4, ..., 1.2,');end
  if(ncase==4)title('ncase = 4, u(z,t), t = 0.2, 0.4, ..., 1.2,');end
  if(ncase==5)title('ncase = 5, u(z,t), t = 0.2, 0.4, ..., 1.2,');end  
%
% Next case
  end
  
