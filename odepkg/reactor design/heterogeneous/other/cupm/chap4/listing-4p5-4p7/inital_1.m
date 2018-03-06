  function u=inital_1(t0)
%
% Function inital_1 is called by the main program to define the
% initial conditions of the four-section retinal O2 transport model
%
% Model parameters
  global   nir     nor     nfl     ncc...
         zl_ir   zl_or   zl_fl   zl_cc...
         zg_ir   zg_or   zg_fl   zg_cc...
          D1ir    D1or    D1fl    D1cc...
               k1ir_or k1or_fl k1fl_cc...
          k1ir    k1or    k1fl    k1cc...
         u1ort    k2or                ...
         pir_s   pcc_s   ncall   ncase      nt
%
% Initial conditions
%
%    Inner retina
     for i=1:nir
      u1ir(i)=pir_s;
      u(i)=u1ir(i);
     end
%
%    Outer retina
     for i=1:nor
       u1or(i)=pir_s;
       u2or(i)=1.0e+07;
       u(i+nir)=u1or(i);
       u(nt+i) =u2or(i);
     end
%
%    Fluid layer
     for i=1:nfl
       u1fl(i)=pir_s;
       u(i+nir+nor)=u1fl(i);
     end
%
%    Choriocapillaris
     for i=1:ncc
       u1cc(i)=pir_s;
       u(i+nir+nor+nfl)=u1cc(i);
     end
% 
% Initialize calls to pde_1
  ncall=0;
 
