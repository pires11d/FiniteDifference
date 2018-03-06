%
% Four-section retinal O2 transport model equation
%
% Clear previous files
  clear all
  clc
%
% Model parameters shared with other routines
  global   nir     nor     nfl     ncc...
         zl_ir   zl_or   zl_fl   zl_cc...
         zg_ir   zg_or   zg_fl   zg_cc...
           Dir     Dor     Dfl     Dcc...
                kir_or  kor_fl  kfl_cc...
           kir     kor     kfl     kcc...           
         pir_s   pcc_s   ncall   ncase 
%
% Select case
%
%   ncase = 1: metabolism rates = 0
%
%   ncase = 2: metabolism rates ne 0
  ncase=2;       
%
% Number of grid points
  nir=11; nor=11; nfl=11; ncc=11;
%
% Lengths of four sections (microns)
  zl_ir=200; zl_or=200; zl_fl=200; zl_cc=200;
%
% Spatial grids;
  zg_ir=[0:zl_ir/10:zl_ir]'; zg_or=[0:zl_or/10:zl_or]';
  zg_fl=[0:zl_fl/10:zl_fl]'; zg_cc=[0:zl_cc/10:zl_cc]';
%
% Diffusivities (microns^2/s)
  Dir=1.0e+04; Dor=1.0e+04; Dfl=1.0e+04; Dcc=1.0e+04;
% Dir=0.1e+04; Dor=0.1e+04; Dfl=0.1e+04; Dcc=0.1e+04;
% Dir=0.5e+04; Dor=0.5e+04; Dfl=0.5e+04; Dcc=0.5e+04;
%
% Interface coefficients
  kir_or=1; kor_fl=1; kfl_cc=1;
%
% Metabolism rates
  if(ncase==1) kir=   0; kor=   0; kfl=   0; kcc=   0; end
  if(ncase==2) kir= 0.1; kor= 0.1; kfl= 0.1; kcc= 0.1; end
%
% Normal breathing O2 concentration (mm hg)
  pir_s=20;
  pcc_s=100;
%
% Independent variable t
  t0=0.0; tf=3.0e+01; tout=[t0:tf/6:tf]'; 
  nout=7;
%  
% Initial condition 
  u0=inital_1(t0);
%
% ODE integration
  mf=2;
  reltol=1.0e-06; abstol=1.0e-05;
  options=odeset('RelTol',reltol,'AbsTol',abstol);
%
% Explicit (nonstiff) integration
  if(mf==1)[t,u]=ode45(@pde_1,tout,u0,options); end
%
% Implicit (sparse stiff) integration
  if(mf==2)
    S=jpattern_num_1;
%   pause
    options=odeset(options,'JPattern',S)
    [t,u]=ode15s(@pde_1,tout,u0,options);   
  end 
%
% One vector to four vectors
  for it=1:nout
    for i=1:nir uir(it,i)=u(it,i);             end
    for i=1:nor uor(it,i)=u(it,i+nir);         end
    for i=1:nfl ufl(it,i)=u(it,i+nir+nor);     end
    for i=1:ncc ucc(it,i)=u(it,i+nir+nor+nfl); end
    if(it>=2) 
        uir(it,1)  =pir_s;
        ucc(it,ncc)=pcc_s; 
    end
  end
%
%   Display tabulated numerical solution
    for it=1:nout
      fprintf('\n    t     zg_ir    uir(z,t)\n')
      for i=1:nir
        fprintf('%5.2f%10.0f%12.2f\n',t(it),zg_ir(i),uir(it,i))
      end  
      fprintf('\n    t     zg_or    uor(z,t)\n')
      for i=1:nor
        fprintf('%5.2f%10.0f%12.2f\n',t(it),zg_or(i),uor(it,i))
      end  
      fprintf('\n    t     zg_fl    ufl(z,t)\n')
      for i=1:nfl
        fprintf('%5.2f%10.0f%12.2f\n',t(it),zg_fl(i),ufl(it,i))
      end  
      fprintf('\n    t     zg_cc    ucc(z,t)\n')
      for i=1:ncc
        fprintf('%5.2f%10.0f%12.2f\n',t(it),zg_cc(i),ucc(it,i))
      end  
%
%   Next t
    end   
    fprintf('\n    ncall = %4d\n',ncall);
%
%   Four plots for uir, uor, ufl, ucc
    figure(2)
    subplot(2,2,1)
    plot(zg_ir,uir,'-')
    axis([0 200 0 100])
    xlabel('zg-ir, \mum'); ylabel('uir(z,t), mm hg O_2');
    title('Inner, t=0,5,...,30 s, D=1\times10^{4}\mum^2/s');  
    subplot(2,2,2)
    plot(zg_or,uor,'-')
    axis([0 200 0 100])
    xlabel('zg-or, \mum'); ylabel('uor(z,t), mm hg O_2'); 
    title('Outer, t=0,10,...,30 s, D=1\times10^{4}\mum^2/s');  
    subplot(2,2,3)
    plot(zg_fl,ufl,'-')
    axis([0 200 0 100])
    xlabel('zg-fl, \mum'); ylabel('ufl(z,t), mm hg O_2');
    title('Fluid, t=0,10,...,30 s, D=1\times10^{4}\mum^2/s');
    subplot(2,2,4)
    plot(zg_cc,ucc,'-')
    axis([0 200 0 100])
    xlabel('zg-cc, \mum'); ylabel('ucc(z,t), mm hg O_2');
    title('Choroid, t=0,10,...,30 s, D=1\times10^{4}\mum^2/s');  
    

    
    
