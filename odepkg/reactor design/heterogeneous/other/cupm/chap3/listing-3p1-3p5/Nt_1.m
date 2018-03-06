  function u0=Nt_1(i)
%
% Function u0 computes the IC for the Nt PDE ATG model 
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
  u0=1.0e+05*(1-tanhr)/2+1.0e+03*(1+tanhr)/2;

 