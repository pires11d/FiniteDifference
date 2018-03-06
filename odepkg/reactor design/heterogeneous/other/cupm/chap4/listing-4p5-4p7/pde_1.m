  function ut=pde_1(t,u)
%
% Function pde_1 defines the ODEs in the MOL solution of the
% four-section retinal O2 transport model equations
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
% One vector to five vectors
%
%   u1
    for i=1:nir u1ir(i)=u(i);                end
    for i=1:nor u1or(i)=u(i+nir);            end
    for i=1:nfl u1fl(i)=u(i+nir+nor);        end
    for i=1:ncc u1cc(i)=u(i+nir+nor+nfl);    end
%
%   u2
    for i=1:nor u2or(i)=u(nt+i);             end    
%
% First order spatial derivatives, u1
  u1ir_z=dss004(0,zl_ir,nir,u1ir); 
  u1or_z=dss004(0,zl_or,nor,u1or); 
  u1fl_z=dss004(0,zl_fl,nfl,u1fl); 
  u1cc_z=dss004(0,zl_cc,ncc,u1cc);   
%
% BCs, u1
%
%   Inner retina
    u1ir(1)=pir_s;     
    u1ir(nir)=k1ir_or*u1or(1);
%
%   Outer retina
    u1or_z(1)=(D1ir/D1or)*u1ir_z(nir);
    u1or(nor)=k1or_fl*u1fl(1); 
%
%   Fluid layer
    u1fl_z(1)=(D1or/D1fl)*u1or_z(nor);
    u1fl(nfl)=k1fl_cc*u1cc(1); 
%
%   Choriocapillaris
    u1cc_z(1)=(D1fl/D1cc)*u1fl_z(nfl);
    if(ncase==1)u1cc(ncc)=pcc_s;    end 
    if(ncase==2)u1cc(ncc)=0.1*pcc_s;end     
%
% Second order spatial derivatives, u1
  nl1=1; nu1=1;
  u1ir_zz=dss044(0,zl_ir,nir,u1ir,u1ir_z,nl1,nu1);
  nl1=2; nu1=1;  
  u1or_zz=dss044(0,zl_or,nor,u1or,u1or_z,nl1,nu1);
  u1fl_zz=dss044(0,zl_fl,nfl,u1fl,u1fl_z,nl1,nu1);
  u1cc_zz=dss044(0,zl_cc,ncc,u1cc,u1cc_z,nl1,nu1);  
%
% PDEs
%
%   u1
    u1ir_t=D1ir*u1ir_zz-k1ir*u1ir;
    u1or_t=D1or*u1or_zz-k1or*u1or;
    u1fl_t=D1fl*u1fl_zz-k1fl*u1fl;  
    u1cc_t=D1cc*u1cc_zz-k1cc*u1cc; 
    u1ir_t(1)  =0;
    u1cc_t(ncc)=0;
%
%   u2
    for i=1:nor
      if(u1or(i)>=u1ort)
        u2or_t(i)=0;
      else
        u2or_t(i)=-k2or*(u1ort-u1or(i));
      end  
    end  
%
% Five vectors to one vector
%
%   u1
    for i=1:nir ut(i)               =u1ir_t(i); end
    for i=1:nor ut(i+nir)           =u1or_t(i); end
    for i=1:nfl ut(i+nir+nor)       =u1fl_t(i); end  
    for i=1:ncc ut(i+nir+nor+nfl)   =u1cc_t(i); end  
%
%   u2
    for i=1:nor ut(nt+i)            =u2or_t(i); end    
  ut=ut';
%
% Increment calls to pde_1
  ncall=ncall+1;
