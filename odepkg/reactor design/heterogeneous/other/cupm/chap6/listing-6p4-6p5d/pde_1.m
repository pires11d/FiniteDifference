  function ut=pde_1(t,u)
%
% Global area
  global      h     cm  alpha   beta...
            u10    Du1     la       ...
            u20    Du2     li       ...
             nr      r  ncase  ncall
%
% One vector to two vectors
  for i=1:nr
    u1(i)=u(i);
    u2(i)=u(i+nr);
  end        
%
% BCs at r = 1
  u1(nr)=u10; u2(nr)=u20; nu=1;
%  
% ur 
  u1r=dss004(0,1,nr,u1);
  u2r=dss004(0,1,nr,u2);  
%
% BCs ar r = 0
  u1r(1)=0; u2r(1)=0; nl=2;  
%
% urr
  u1rr=dss044(0,1,nr,u1,u1r,nl,nu);
  u2rr=dss044(0,1,nr,u2,u2r,nl,nu);  
%
% ut
  for i=1:nr-1
    if(ncase==1) lam=la; sc=sa(u2(i)); g=ga(u1(i)); end
    if(ncase==2) lam=li; sc=si(u2(i)); g=gi(u1(i)); end
    if(i==1) 
      u1t(i)=Du1*2*u1rr(i)+sc*u1(i)*(2-u1(i))-u1(i);
      u2t(i)=Du2*2*u2rr(i)+lam*g-lam*u2(i); 
    else  
    u1t(i)=Du1*(u1rr(i)+(1/r(i))*u1r(i))+sc*u1(i)*(2-u1(i))-u1(i);
    u2t(i)=Du2*(u2rr(i)+(1/r(i))*u2r(i))+lam*g-lam*u2(i);    
    end
  end  
  u1t(nr)=0; u2t(nr)=0;
%
% Two vectors to one vector
  for i=1:nr
    ut(i)   =u1t(i);
    ut(i+nr)=u2t(i);
  end  
%
% Transpose and counter
  ut=ut';
  ncall=ncall+1;