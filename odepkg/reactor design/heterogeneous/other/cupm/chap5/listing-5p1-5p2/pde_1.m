  function ut=pde_1(u,t)
%
% Function pde_1 computes the t derivative vector of the u vector 
%
  global zL q1 q2 v1 v2 V1L V1R kM1 kM2 u1L u2zL u1 u2 ifd n ncase ncall
%
% One vector to two PDEs, one ODE
  for i=1:n
    u1(i)=u(i);
    u2(i)=u(i+n);
  end 
  u1R=u(2*n+1);
%
% Boundary condition
  u2(n)=u2zL;
%
% First order spatial derivative
%
%   ifd = 1: Two point upwind finite difference (2pu)
    if(ifd==1) [u1z]=dss012(0.0,zL,n,u1,v1); end
    if(ifd==1) [u2z]=dss012(0.0,zL,n,u2,v2); end    
%      
%   ifd = 2: Three point center finite difference (3pc)
    if(ifd==2) [u1z]=dss002(0.0,zL,n,u1); end
    if(ifd==2) [u2z]=dss002(0.0,zL,n,u2); end    
%      
%   ifd = 3: Five point biased upwind approximation (5pbu)      
    if(ifd==3) [u1z]=dss020(0.0,zL,n,u1,v1); end  
    if(ifd==3) [u2z]=dss020(0.0,zL,n,u2,v2); end      
%      
%   ifd = 4: van Leer flux limiter  
    if(ifd==4) [u1z]=vanl(0.0,zL,n,u1,v1); end 
    if(ifd==4) [u2z]=vanl(0.0,zL,n,u2,v2); end 
%      
%   ifd = 5: Superbee flux limiter     
    if(ifd==5) [u1z]=super(0.0,zL,n,u1,v1); end
    if(ifd==5) [u2z]=super(0.0,zL,n,u2,v2); end    
%      
%   ifd = 6: Smart flux limiter      
    if(ifd==6) [u1z]=smart(0.0,zL,n,u1,v1); end 
    if(ifd==6) [u2z]=smart(0.0,zL,n,u2,v2); end 
%
% Temporal derivatives
%
%    u1t
     u1t(1)=(1/V1L)*q1*(u1L-u1(1));
     for i=2:n
       u1t(i)=-v1*u1z(i)+kM1*(u2(i)-u1(i));
     end
     u1Rt=(1/V1R)*q1*(u1(n)-u1R);
%
%    u2t
     u2t(n)=0.0;
     for i=1:n-1
       u2t(i)=-v2*u2z(i)+kM2*(u1(i)-u2(i));
     end 
%
% Two PDEs, one ODE to one vector
  for i=1:n
      ut(i)=u1t(i);
    ut(i+n)=u2t(i);
  end  
  ut(2*n+1)=u1Rt;
%
% Increment calls to pde_1
  ncall=ncall+1;
