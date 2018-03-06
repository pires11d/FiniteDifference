  function ut=pde_1(t,u)
%
% Global area
  global     nr     r0      r     u0...
              D      p     sc  ncall
%
% BC at r = 0
  u(1)=u0; nl=1;
%  
% ur 
  ur=dss004(0,r0,nr,u);
%
% BC ar r = r0
  ur(nr)=0; nu=2;  
%
% urr
  urr=dss044(0,r0,nr,u,ur,nl,nu);
%
% ut
  for i=2:nr
    ut(i)=D*((1-u(i)/u0)^p*urr(i)+p*(1-u(i)/u0)^(p-1)*(-1/u0)*ur(i)^2)+...
          sc*u(i)*(1-u(i)/u0);
  end
  ut(1)=0;
%
% Transpose and counter
  ut=ut';
  ncall=ncall+1;