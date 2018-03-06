  function ut=pde_1(u,t)
%
% Function pde_1 computes the t derivative vector of the u vector 
%
  global zL D ifd n ncase ncall
%
% Second order spatial derivative
  if(ncase==1)
     [uz]=dss002(0,zL,n,u);
     uz(1)=u(1); uz(n)=-u(n);
    [uzz]=dss002(0,zL,n,uz);
  end 
  if(ncase==2)
     [uz]=dss004(0,zL,n,u);
     uz(1)=u(1); uz(n)=-u(n);     
    [uzz]=dss004(0,zL,n,uz);
  end  
  if(ncase==3)
     [uz]=dss006(0,zL,n,u);
     uz(1)=u(1); uz(n)=-u(n);     
    [uzz]=dss006(0,zL,n,uz);
  end 
  if(ncase==4)
    nl=2; nu=2; uz(1)=u(1); uz(n)=-u(n);  
    [uzz]=dss042(0,zL,n,u,uz,nl,nu);
  end 
  if(ncase==5)
    nl=2; nu=2; uz(1)=u(1); uz(n)=-u(n);       
    [uzz]=dss044(0,zL,n,u,uz,nl,nu);
  end 
  if(ncase==6)
    nl=2; nu=2; uz(1)=u(1); uz(n)=-u(n);       
    [uzz]=dss046(0,zL,n,u,uz,nl,nu);
  end 
%
% Temporal derivative
  ut=D*uzz;
% for i=1:n
%   ut(i)=D*uzz(i);
% end
%
% Increment calls to pde_1
  ncall=ncall+1;
