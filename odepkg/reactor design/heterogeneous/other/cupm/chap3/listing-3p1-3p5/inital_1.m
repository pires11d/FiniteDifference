  function u=inital_1(t0)
%
% Function inital_1 is called by the main program to define the
% initial conditions of the 3-PDE ATG model
%
% Model parameters
  global    Nn      Nt      Ch...
            rl      ru       r       n...
           rn1     rn2      Kn...
           rt1      Dt      Kt...
           rh1     rh2      Dh...
         ncall   ncase    
%
% Initial conditions
  for i=1:n
    u(i)    =Nn_1(i);
    u(i+n)  =Nt_1(i);
    u(i+2*n)=Ch_1(i);
  end
% 
% Initialize calls to inital_1
  ncall=0;
 
