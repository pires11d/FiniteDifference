  clc
  clear all
% 
%   One dimensional diffusion
%
%   The PDE that models the diffusion mass transfer is
%
%      ut = D*uzz                                              (1)
%
%   with the initial and Neumann boundary conditions
%
%      u(z,0) = u0(z), uz(0,t) = 0, uz(zL,t) = 0         (2)(3)(4)
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
%   The profiles of u(z,t) in z are plotted for selected t.
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
  D=1; zL=1; n=21;
%
% Step through cases
  for ncase=1:6
%
% Parameters for fourth order Runge Kutta integration
  nsteps=100;
  h=0.0005;
%
% Initial condition
  for i=1:n
    z(i)=(i-1)/(n-1)*zL;  
    u(i)=cos(pi*z(i)/zL);
  end
  t=0;
%
% Write ncase, h
  fprintf('\n\n ncase = %5d    h = %10.3e\n\n',ncase,h);
%
% Analytical solution at t = 0
  for i=1:n
    ua(i)=cos(pi*z(i)/zL);
%
%   Store solutions for plotting at t = 0
     uplot(1,i)=u(i);
    uaplot(1,i)=ua(i);
     zplot(i)=z(i);
  end
%
% nout output points
  nout=7;
  ncall=0;
  for iout=2:nout
%
%   Fourth order Runge Kutta integration
    u0=u; t0=t;
    [u,t]=rk4(u0,t0,h,nsteps);      
%
%   Analytical solution 
    for i=1:n
      ua(i)=exp(-(pi/zL)^2*t)*cos(pi*z(i)/zL);
%
%     Store solutions for plotting
       uplot(iout,i)=u(i);
      uaplot(iout,i)=ua(i);
%
%   Next point in z
    end
%
% Next output
  end    
%
% Plot z profiles of solution
  figure(ncase);
  plot(zplot,uplot(1,:),'-o',zplot,uaplot(1,:),'-');     
  axis([0 1 -1 1]);
  ylabel('u(z,t),ua(z,t)');xlabel('t');
  if(ncase==1)title('ncase = 1; t = 0, 0.05,..., 0.3; num - o; anal - line');end
  if(ncase==2)title('ncase = 2; t = 0, 0.05,..., 0.3; num - o; anal - line');end
  if(ncase==3)title('ncase = 3; t = 0, 0.05,..., 0.3; num - o; anal - line');end
  if(ncase==4)title('ncase = 4; t = 0, 0.05,..., 0.3; num - o; anal - line');end
  if(ncase==5)title('ncase = 5; t = 0, 0.05,..., 0.3; num - o; anal - line');end
  if(ncase==6)title('ncase = 6; t = 0, 0.05,..., 0.3; num - o; anal - line');end
  hold on
  for iout=2:nout
    plot(zplot,uplot(iout,:),'-o',zplot,uaplot(iout,:),'-');        
    hold on
  end  
%
% Next case
  end

