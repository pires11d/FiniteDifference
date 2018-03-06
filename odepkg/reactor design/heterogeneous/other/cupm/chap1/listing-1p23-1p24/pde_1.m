  function ut=pde_1(u,t)
%
% Function pde_1 computes the t derivative vector of the u vector 
%
  global a b D v ue kr zl n ncase ncall
%
% First, second order spatial derivatives
  nl=1; nu=2;
  if(ncase==1)
     [uz]=dss012(0,zl,n,u,v);
     u(1)=ue-D*uz(1); 
     [uzz]=dss042(0,zl,n,u,uz,nl,nu);
  end 
  if(ncase==2)
     [uz]=dss002(0,zl,n,u);
     u(1)=ue-D*uz(1); 
     [uzz]=dss042(0,zl,n,u,uz,nl,nu);
  end 
  if(ncase==3)
     [uz]=vanl(0,zl,n,u,v);
     u(1)=ue-D*uz(1); 
     [uzz]=dss042(0,zl,n,u,uz,nl,nu);
  end 
  if(ncase==4)
     [uz]=super(0,zl,n,u,v);      
     u(1)=ue-D*uz(1); 
     [uzz]=dss044(0,zl,n,u,uz,nl,nu);     
  end  
  if(ncase==5)
     [uz]=smart(0,zl,n,u,v);
     u(1)=ue-D*uz(1); 
     [uzz]=dss046(0,zl,n,u,uz,nl,nu);     
  end 
%
% Temporal derivative
  for i=2:n-1
    ut(i)=D*uzz(i)-v*uz(i)-kr*u(i);
  end
  ut(1)=0; ut(n)=-v*uz(n);
%
% Increment calls to pde_1
  ncall=ncall+1;
