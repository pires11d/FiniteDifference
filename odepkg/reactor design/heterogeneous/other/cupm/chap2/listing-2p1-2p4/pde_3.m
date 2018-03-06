  function ut=pde_3(t,u)
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
% BCs
  c(n)=cbulk;
  cz(1)=(1/D)*(kf*c(1)*(cbsat-cb)-kr*cb);
  nu=1; nl=2;
%
% czz
  czz=dss044(zl,zu,n,c,cz,nl,nu);
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
% Increment calls to pde_3
  ncall=ncall+1;

  
