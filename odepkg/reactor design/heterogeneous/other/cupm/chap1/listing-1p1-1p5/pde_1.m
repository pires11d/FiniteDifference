  function ut=pde_1(u,t)
%
% Function pde_1 computes the t derivative vector of the u vector 
%
  global z dz zL v n ncase ncall
%
% Boundary condition
  if(ncase==1) u(1)=1; end
  if(ncase==2) u(n)=1; end
%
% Temporal derivative
  if(ncase==1)
    for i=2:n
      ut(i)=-v*(u(i)-u(i-1))/dz;
    end 
    ut(1)=0;
  end
  if(ncase==2)
    for i=1:n-1
      ut(i)=-v*(u(i+1)-u(i))/dz;
    end 
    ut(n)=0;
  end
%
% Increment calls to pde_1
  ncall=ncall+1;
