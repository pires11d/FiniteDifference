  function u0=Nn_1(i)
%
% Function u0 computes the IC for the Nn PDE ATG model 
%
  global    Nn      Nt      Ch...
            rl      ru       r       n...
           rn1     rn2      Kn...
           rt1      Dt      Kt...
           rh1     rh2      Dh...
         ncall   ncase
%
% IC function
  rs=50;
  num=exp(rs*(r(i)-r(21)))-exp(-rs*(r(i)-r(21)));
  den=exp(rs*(r(i)-r(21)))+exp(-rs*(r(i)-r(21)));
  tanhr=num/den;
  u0=5.0e+07*(1-tanhr)/2+1.0e+08*(1+tanhr)/2;


 