  function ut=pde_1(u,t)
%
% Function pde_1 computes the t derivative vector of the u vector 
%
  global zL D ifd n ncase ncall
%
% Boundary conditions at z = 0, z=zL
  u(1)=0; u(n)=0;
%
% Second order spatial derivative
  if(ncase==1)
     [uz]=dss002(0,zL,n,u);
    [uzz]=dss002(0,zL,n,uz);
  end 
  if(ncase==2)
     [uz]=dss004(0,zL,n,u);
    [uzz]=dss004(0,zL,n,uz);
  end  
  if(ncase==3)
     [uz]=dss006(0,zL,n,u);
    [uzz]=dss006(0,zL,n,uz);
  end 
  if(ncase==4)
    nl=1; nu=1; uz(1)=0;  
    [uzz]=dss042(0,zL,n,u,uz,nl,nu);
  end 
  if(ncase==5)
    nl=1; nu=1; uz(1)=0; 
    [uzz]=dss044(0,zL,n,u,uz,nl,nu);
  end 
  if(ncase==6)
    nl=1; nu=1; uz(1)=0;  
    [uzz]=dss046(0,zL,n,u,uz,nl,nu);
  end 
%
% Temporal derivative
  ut(1)=0.0; ut(n)=0;
  for i=2:n-1
    ut(i)=D*uzz(i);
  end
%
% Increment calls to pde_1
  ncall=ncall+1;
