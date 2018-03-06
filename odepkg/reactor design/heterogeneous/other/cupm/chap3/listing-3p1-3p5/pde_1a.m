  function ut=pde_1(t,u)
%
% Function pde_1 defines the ODEs in the MOL solution of three
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
    if(i==1)
      Nnt(i)=rn1*Nn(i)*(1-Nn(i)/Kn)-rn2*Ch(i)*Nn(i);
      Ntt(i)=rt1*Nt(i)*(1-Nt(i)/Kt)+3*Dt*Ntrr(i);
      Cht(i)=rh1*Nt(i)-rh2*Ch(i)+3*Dh*Chrr(i);
    else 
      Nnt(i)=rn1*Nn(i)*(1-Nn(i)/Kn)-rn2*Ch(i)*Nn(i);
      Ntt(i)=rt1*Nt(i)*(1-Nt(i)/Kt)+Dt*(Ntrr(i)+2/r(i)*Ntr(i));
      Cht(i)=rh1*Nt(i)-rh2*Ch(i)+Dh*(Chrr(i)+2/r(i)*Chr(i));
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
% Increment calls to pde_1
  ncall=ncall+1;
