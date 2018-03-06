  function ut=pde_1(u,t)
%
% Function pde_1 computes the t derivative vector of the u vector 
%
  global z dz zL v ifd n ncase ncall
%
% Boundary condition
  if(ncase==1) u(1)=1; end
  if(ncase==2) u(n)=1; end
%
% uz; approximation selected by ifd
%
%   ifd = 1: Two point upwind finite difference (2pu)
    if(ifd==1) [uz]=dss012(0,zL,n,u,v); end
%      
%   ifd = 2: Three point center finite difference (3pc)
    if(ifd==2) [uz]=dss002(0,zL,n,u); end
%      
%   ifd = 3: Five point biased upwind approximation (5pbu)      
    if(ifd==3) [uz]=dss020(0,zL,n,u,v); end  
%      
%   ifd = 4: van Leer flux limiter  
    if(ifd==4) [uz]=vanl(0,zL,n,u,v); end 
%      
%   ifd = 5: Superbee flux limiter     
    if(ifd==5) [uz]=super(0,zL,n,u,v); end  
%      
%   ifd = 6: Smart flux limiter      
    if(ifd==6) [uz]=smart(0,zL,n,u,v); end         
%
% ut (PDE)
  ut=-v*uz;
  if(ncase==1) ut(1)=0; end
  if(ncase==2) ut(n)=0; end
%
% Increment calls to pde_1
  ncall=ncall+1;
