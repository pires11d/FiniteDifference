  function ut=pde_1(t,u)
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
  cf=c(2)-(2*dz/D)*(kf*c(1)*(cbsat-cb)-kr*cb);
  c(n)=cbulk;
%
% PDE
  for i=1:n
    if(i==1)     ct(1)=D*(c(2)-2*c(1)+cf)/dz2; 
    elseif(i==n) ct(n)=0;
    else         ct(i)=D*(c(i+1)-2*c(i)+c(i-1))/dz2;
    end    
  end 
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
% Increment calls to pde_1
  ncall=ncall+1;

  
