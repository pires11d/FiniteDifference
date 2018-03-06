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
          D2ir    D2or    D2fl    D2cc...          
               k1ir_or k1or_fl k1fl_cc...
               k2ir_or k2or_fl k2fl_cc...               
          k1ir    k1or    k1fl    k1cc...
          k2ir    k2or    k2fl    k2cc...
         u1irt   u1ort   u1flt   u1cct... 
         pir_s   pcc_s   ncall   ncase      nt
%
% Initial conditions
%
%    Inner retina
     for i=1:nir
      u1ir(i)=pir_s;
      u2ir(i)=0;
        u(i)   =u1ir(i);
        u(nt+i)=u2ir(i);
     end
%
%    Outer retina
     for i=1:nor
       u1or(i)=pir_s;
       u2or(i)=0;
       u(i+nir)   =u1or(i);
       u(nt+i+nir)=u2or(i);
     end
%
%    Fluid layer
     for i=1:nfl
       u1fl(i)=pir_s;
       u2fl(i)=0;
       u(i+nir+nor)   =u1fl(i);
       u(nt+i+nir+nor)=u2fl(i);
     end
%
%    Choriocapillaris
     for i=1:ncc
       u1cc(i)=pir_s;
       u2cc(i)=0;
       u(i+nir+nor+nfl)   =u1cc(i);
       u(nt+i+nir+nor+nfl)=u2cc(i);
     end
% 
% Initialize calls to pde_1
  ncall=0;
 
