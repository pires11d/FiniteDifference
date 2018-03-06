  function u=inital_1(t0)
%
% Function inital_1 is called by the main program to define the
% initial conditions of the four-section retinal O2 transport model
%
% Model parameters
  global   nir     nor     nfl     ncc...
         zl_ir   zl_or   zl_fl   zl_cc...
         zg_ir   zg_or   zg_fl   zg_cc...
           Dir     Dor     Dfl     Dcc...
                kir_or  kor_fl  kfl_cc...
           kir     kor     kfl     kcc...
         pir_s   pcc_s   ncall   ncase 
%
% Initial conditions
%
%    Inner retina
     for i=1:nir
      uir(i)=0;
        u(i)=uir(i);
     end
%
%    Outer retina
     for i=1:nor
       uor(i)=0;
       u(i+nir)=uor(i);
     end
%
%    Fluid layer
     for i=1:nfl
       ufl(i)=0;
       u(i+nir+nor)=ufl(i);
     end
%
%    Choriocapillaris
     for i=1:ncc
       ucc(i)=0;
       u(i+nir+nor+nfl)=ucc(i);
     end
% 
% Initialize calls to pde_1
  ncall=0;
 
