  function ut=pde_2(t,u)
%
% Problem parameters
  global zl zu z dz dz2 D kf cbsat kr n cbulk ndss ncall
%
% ODE and PDE
  for i=1:n
    c(i)=u(i);
  end
  cb=u(n+1);
%
% BC
  c(n)=cbulk;
%
% cz
  cz=dss004(zl,zu,n,c);
%
% BC
  cz(1)=(1/D)*(kf*c(1)*(cbsat-cb)-kr*cb);
%
% czz
  czz=dss004(zl,zu,n,cz);
%
% PDE
  ct=D*czz;
  ct(n)=0;  
%
% ODE
  cbt=kf*c(1)*(cbsat-cb)-kr*cb;
%
% Derivative vector
  for i=1:n
    ut(i)=ct(i);
  end
  ut(n+1)=cbt;
%
% Transpose for ODE integrator
  ut=ut';
%
% Increment calls to pde_2
  ncall=ncall+1;

  
