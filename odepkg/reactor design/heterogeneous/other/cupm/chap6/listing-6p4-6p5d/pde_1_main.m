%
% Wound healing
%
% Clear previous files
  clear all
  clc
%
% Global area
  global      h     cm  alpha   beta...
            u10    Du1     la       ...
            u20    Du2     li       ...
             nr      r  ncase  ncall
%
% Case
%
%   ncase = 1: Activator
%
%   ncase = 2: Inhibitor
  ncase=2;
%
% Model parameters
  la=30; li=5; h=10; cm=40; alpha=0.1; 
  beta=(1+cm^2-2*h*cm)/(1-cm)^2;
  if(ncase==1) Du1=5.0e-04; Du2=0.45; end  
  if(ncase==2) Du1=1.0e-04; Du2=0.85; end
  fprintf('\n ncase = %2d\n',ncase);
  fprintf('\n  la = %5.1f     li = %5.1f      h = %5.1f\n',la,li,h);   
  fprintf('\n  cm = %5.1f  alpha = %5.3f   beta = %5.3f\n',cm,alpha,beta);
  fprintf('\n Du1 = %8.3e    Du2 = %8.3e\n',Du1,Du2);
%
% Spatial grid
  nr=101; dr=1/(nr-1); 
  for j=1:nr
    r(j)=(j-1)*dr;
  end
  u10=1; u20=1;
%
% Independent variable for ODE integration
  if(ncase==1) tf=25; end
  if(ncase==2) tf=50; end
  tout=[0:tf/20:tf]'; nout=21;  
  ncall=0;
%
% Initial condition
  for i=1:nr u0(i)=0; u0(i+nr)=0; end
%
% ODE integration 
  reltol=1.0e-04; abstol=1.0e-04;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
  [t,u]=ode15s(@pde_1,tout,u0,options);   
%
% Display a heading and t at output
  fprintf('\n  nr = %2d\n',nr);
  for it=1:nout
    for i=1:nr
      u1(it,i)=u(it,i);
      u2(it,i)=u(it,nr+i);
    end
    u1(it,nr)=u10; u2(it,nr)=u20;    
    fprintf('\n t = %6.2e',t(it));
%
%   Detailed solution
%   for i=1:nr
%     fprintf('\n i = %3d  %5.3f  %5.3f  %5.3f  ',...
%             i,r(i),u1(it,i),u2(it,i));
%   end      
  end
  fprintf('\n\n ncall = %5d\n',ncall);
%
% Plot the solution as u1(r,t), u2(r,t) vs r with t as a parameter  
  figure(1)
  plot(r,u1(2:nout,:),'-');
  xlabel('r'); ylabel('u1(r,t)');
  title('u1(r,t)vs r')
  figure(2)
  plot(r,u2(2:nout,:),'-');
  xlabel('r'); ylabel('u2(r,t)');
  title('u2(r,t) vs r') 
  
  
