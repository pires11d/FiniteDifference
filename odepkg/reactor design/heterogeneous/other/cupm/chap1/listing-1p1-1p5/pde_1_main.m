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
  global z dz zL v n ncase ncall
%
% Grid (in z)
  zL=1; n=51; dz=0.02;
  z=[0:dz:zL];
%
% Level of output
%
%   Detailed output - ip = 1
%
%   Brief (IC) output - ip = 2
%
  ip=2;
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
  nsteps=20;
  h=0.001;
%
% Initial condition
  for i=1:n
    u(i)=0;
  end
  t=0;
%
% Write ncase, h, v
  fprintf('\n\n ncase = %5d    h = %10.3e   v = %4.2f\n\n',ncase,h,v);
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
% Display numerical, analytical solutions at t = 0
  if(t <  zL/abs(v)) ua=0;   end
  if(t >  zL/abs(v)) ua=1;   end
  if(t == zL/abs(v)) ua=0.5; end 
  if(ncase==1)
    diff=u(n)-ua;
    fprintf('%5.2f%7.2f%10.2f%10.3f%10.3f%10.4f\n',...
            t,zL,t-zL/abs(v),u(n),ua,diff);
  end 
  if(ncase==2)
    diff=u(1)-ua;
    fprintf('%5.2f%7.2f%10.2f%10.3f%10.3f%10.4f\n',...
            t,zL,t-zL/abs(v),u(1),ua,diff);
  end  
%
% Store solution for plotting
  if(ncase==1)
     uplot(1,1)=u(n);
    uaplot(1,1)=ua;
     tplot(1)=t;
  end 
  if(ncase==2)
     uplot(2,1)=u(1);
    uaplot(2,1)=ua;
  end   
%
% nout output points
  nout=101;
  ncall=0;
  for iout=2:nout
%
%   Euler integration
    u0=u; t0=t;
    [u,t]=euler(u0,t0,h,nsteps);
%
%   Numerical, analytical solutions
    if(t <  zL/abs(v)) ua=0;   end
    if(t >  zL/abs(v)) ua=1;   end
    if(t == zL/abs(v)) ua=0.5; end 
    if(ip==1)
      if(ncase==1)
        diff=u(n)-ua;
        fprintf('%5.2f%7.2f%10.2f%10.3f%10.3f%10.4f\n',...
                t,zL,t-zL/abs(v),u(n),ua,diff);
      end 
      if(ncase==2)
        diff=u(1)-ua;
        fprintf('%5.2f%7.2f%10.2f%10.3f%10.3f%10.4f\n',...
                t,zL,t-zL/abs(v),u(1),ua,diff);
      end  
    end
%
%   Store solution for plotting
    if(ncase==1)
       uplot(1,iout)=u(n);
      uaplot(1,iout)=ua;
       tplot(iout)=t;    
    end  
    if(ncase==2)
       uplot(2,iout)=u(1);
      uaplot(2,iout)=ua;
       tplot(iout)=t;    
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
  end

  

