  clc
  clear all
% 
%   One dimensional diffusion
%
%   The PDE that models the diffusion mass transfer is
%
%      ut = D*uzz                                              (1)
%
%   with the initial and Dirichlet boundary conditions
%
%      u(z,0) = u0(z), u(0,t) = 0, u(zL,t) = 0           (2)(3)(4)
%
%   where
%
%      u     dependent variable, e.g., concentration
%
%      t     time
%
%      z     position
%
%      D     diffusivity
%
%      u0(z) initial u
%
%      zL    length
%
%   Of particular interest is the concentration, u(z=zL/2,t).
%
%   The method of lines (MOL) solution of eq. (1) is coded below.
%   The resulting system of ODEs in t, defined on a 21-point grid in z, 
%   is then integrated by the classical fourth order Runge Kutta method.
%
%   The analytical solution to eqs. (1) to (4) is also programmed
%   for comparison with the numerical solution.
%
%   The following code is run for the six special cases:
%
%      Case 1: uzz by three point finite differences (3pc)
%              and stagewise differentiation
%
%      Case 2: uzz by five point finite differences  (5pc)
%              and stagewise differentiation
%
%      Case 3: uzz by seven point finite differences (7pc)
%              and stagewise differentiation
%
%      Case 4: uzz by three point finite differences (3pc)
%              and direct second order differentiation
%
%      Case 5: uzz by five point finite differences  (5pc)
%              and direct second order differentiation
%
%      Case 6: uzz by seven point finite differences (7pc)
%              and direct second order differentiation
%
  global zL D ncase n ncall
%
% Model parameters (defined in the comments above)
  D=1; zL=1; n=21; n2=11;
%
% Step through cases
  for ncase=1:6
%
% Level of output
%
%   Detailed output - ip = 1
%
%   Brief (IC) output - ip = 2
%
  ip=1;
%
% Parameters for fourth order Runge Kutta integration (Note: for 
% nsteps = 10, h = 0.002, the numerical solution is unstable for 
% ncase = 4,5,6)
  nsteps=20;
  h=0.001;
%
% Initial condition
  for i=1:n
    z(i)=(i-1)/(n-1)*zL;  
    u(i)=sin(pi*z(i)/zL);
  end
  t=0;
%
% Write ncase, h
  fprintf('\n\n ncase = %5d    h = %10.3e\n\n',ncase,h);
%
% Write heading
  fprintf('   t     u(zL/2,t)   ua(zL/2,t)      diff\n');
%
% Display numerical, analytical solutions at t = 0, z = zL/2
  ua=sin(pi*z(n2)/zL);
  diff=u(n2)-ua;
  fprintf('%5.2f%12.5f%12.5f%15.4e\n',t,u(n2),ua,diff);
%
% Store solution for plotting
   uplot(ncase,1)=u(n2);
  uaplot(ncase,1)=ua;
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
%   Numerical, analytical solutions
    ua=exp(-(pi/zL)^2*t)*sin(pi*z(n2)/zL);
    diff=u(n2)-ua;
    if(ip==1)
      fprintf('%5.2f%12.5f%12.5f%15.4e\n',t,u(n2),ua,diff);
    end  
%
%   Store solution for plotting
     uplot(ncase,iout)=u(n2);
    uaplot(ncase,iout)=ua;
     tplot(iout)=t;    
%
% Next output
  end
%
% Plots for ncase = 1 to 6
  figure(ncase);
  plot(tplot,uplot(1,:),'-o');
  axis([0 1 0 1]);
  ylabel('u(z=zL/2,t),ua(z=zL/2,t)');xlabel('t');
  if(ncase==1)title('ncase = 1; num - o; anal - line');end
  if(ncase==2)title('ncase = 2; num - o; anal - line');end
  if(ncase==3)title('ncase = 3; num - o; anal - line');end
  if(ncase==4)title('ncase = 4; num - o; anal - line');end
  if(ncase==5)title('ncase = 5; num - o; anal - line');end
  if(ncase==6)title('ncase = 6; num - o; anal - line');end
  hold on
  plot(tplot,uaplot(1,:),'-');  
%
% Next case
  end

