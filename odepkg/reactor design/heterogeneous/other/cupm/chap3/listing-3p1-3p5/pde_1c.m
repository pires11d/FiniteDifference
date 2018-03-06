  function [ut,term]=pde_1c(t,u)
%
% Function pde_1c defines the RHS terms and ODEs in the MOL solution of three
% PDEs for the ATG model
%
  global    Nn      Nt      Ch...
            rl      ru       r       n...
           rn1     rn2      Kn...
           rt1      Dt      Kt...
           rh1     rh2      Dh...
         ncall   ncase
%
% One vector to three vectors
  for i=1:n
    Nn(i)=u(i);
    Nt(i)=u(i+n);
    Ch(i)=u(i+2*n);
  end
%
% First order spatial derivatives
  Nnr=dss004(rl,ru,n,Nn); 
  Ntr=dss004(rl,ru,n,Nt);
  Chr=dss004(rl,ru,n,Ch); 
%
% BCs
  Ntr(1)=0; Ntr(n)=0;
  Chr(1)=0; Chr(n)=0;    
%
% Second order spatial derivatives
  nl=2; nu=2;
  Ntrr=dss044(rl,ru,n,Nt,Ntr,nl,nu);
  Chrr=dss044(rl,ru,n,Ch,Chr,nl,nu);  
%
% PDEs
  for i=1:n
    D=Dt*(1-Nn(i)/Kn);
    if(D<0)D=0;end
    if(i==1)
      term(1,i)=rn1*Nn(i)*(1-Nn(i)/Kn);
      term(2,i)=-rn2*Ch(i)*Nn(i);
      Nnt(i)=term(1,i)+term(2,i);
      term(3,i)=rt1*Nt(i)*(1-Nt(i)/Kt);
      term(4,i)=+3*D*Ntrr(i);
      Ntt(i)=term(3,i)+term(4,i);
      term(6,i)=rh1*Nt(i);
      term(7,i)=-rh2*Ch(i);
      term(8,i)=+3*Dh*Chrr(i);
      Cht(i)=term(6,i)+term(7,i)+term(8,i);
    else 
      term(1,i)=rn1*Nn(i)*(1-Nn(i)/Kn);
      term(2,i)=-rn2*Ch(i)*Nn(i);
      Nnt(i)=term(1,i)+term(2,i);
      term(3,i)=rt1*Nt(i)*(1-Nt(i)/Kt);
      term(4,i)=+D*(Ntrr(i)+2/r(i)*Ntr(i));
      term(5,i)=+(-Dt/Kn)*Nnr(i)*Ntr(i);
      Ntt(i)=term(3,i)+term(4,i)+term(5,i);
      term(6,i)=rh1*Nt(i);
      term(7,i)=-rh2*Ch(i);
      term(8,i)=+Dh*(Chrr(i)+2/r(i)*Chr(i));
      Cht(i)=term(6,i)+term(7,i)+term(8,i);
    end
  end
%
% Three vectors to one vector
  for i=1:n
    ut(i)    =Nnt(i);
    ut(i+n)  =Ntt(i);
    ut(i+2*n)=Cht(i);
  end
  ut=ut';
%
% Increment calls to pde_1c
  ncall=ncall+1;
