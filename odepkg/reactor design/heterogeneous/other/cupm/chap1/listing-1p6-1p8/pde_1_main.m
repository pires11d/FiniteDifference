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
%    t (temporal, initial value) integration: ode45 or ode15s
% 
  global dz zL v n ncase ncall
%
% Grid (in z)
  zL=1; n=51; dz=0.02;
%
% Level of output
%
%   Detailed output - ip = 1
%
%   Brief (IC) output - ip = 2
%
  ip=1;
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
% Write ncase, v
  fprintf('\n\n ncase = %5d  v = %4.2f\n\n',ncase,v);
%
% Write heading
  if(ncase==1)
    fprintf('   t     zL    t-zL/|v|   u(zL,t)   ua(zL,t)    diff\n');
  end  
%
% Write heading
  if(ncase==2)
    fprintf('   t     zL    t-zL/|v|    u(0,t)    ua(0,t)    diff\n');
  end   
%
% Initial condition
  for i=1:n
    u0(i)=0;
  end
  t=0;     
%
% Independent variable for ODE integration
  t0=0;
  tf=2;
  tout=[t0:0.02:tf];
  nout=101;
  ncall=0;
%
% ODE integration
  mf=2;
  reltol=1.0e-06; abstol=1.0e-06;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
%
% Explicit (nonstiff) integration
  if(mf==1)[t,u]=ode45(@pde_1,tout,u0,options); end
%
% Implicit (sparse stiff) integration
  if(mf==2)
    S=jpattern_num;
    options=odeset(options,'JPattern',S)
    [t,u]=ode15s(@pde_1,tout,u0,options);   
  end  
%
% nout output points
  for iout=1:nout
%
%   Numerical, analytical solutions
    if(t(iout) <  zL/abs(v)) ua=0;   end
    if(t(iout) >  zL/abs(v)) ua=1;   end
    if(t(iout) == zL/abs(v)) ua=0.5; end 
    if(ip==1)
      if(ncase==1)
        diff=u(iout,n)-ua;
        fprintf('%5.2f%7.2f%10.2f%10.3f%10.3f%10.4f\n',...
                t(iout),zL,t(iout)-zL/abs(v),u(iout,n),ua,diff);
      end 
      if(ncase==2)
        diff=u(iout,1)-ua;
        fprintf('%5.2f%7.2f%10.2f%10.3f%10.3f%10.4f\n',...
                t(iout),zL,t(iout)-zL/abs(v),u(iout,1),ua,diff);
      end  
    end
%
%   Store solution for plotting
    if(ncase==1)
       uplot(1,iout)=u(iout,n);
      uaplot(1,iout)=ua;
       tplot(iout)=t(iout);    
    end  
    if(ncase==2)
       uplot(2,iout)=u(iout,1);
      uaplot(2,iout)=ua;
       tplot(iout)=t(iout);    
    end   
%
% Next output
  end
%
% Plots for ncase = 1, 2
  if(ncase==1)
    figure(1);
    plot(tplot,uplot(1,:),'-o');
    axis([0 2 0 1]);
    ylabel('u(zL,t),ua(zL,t)');xlabel('t');
    title('ncase = 1; num - o; anal - line');
    hold on
    plot(tplot,uaplot(1,:),'-');
  end  
  if(ncase==2)
    figure(2);
    plot(tplot,uplot(2,:),'-o');
    axis([0 2 0 1]);
    ylabel('u(0,t),ua(0,t)');xlabel('t');
    title('ncase = 2; num - o; anal - line');
    hold on
    plot(tplot,uaplot(2,:),'-');
  end    
%
% Next case
  fprintf('\n ncall = %4d\n',ncall);
  end

  

