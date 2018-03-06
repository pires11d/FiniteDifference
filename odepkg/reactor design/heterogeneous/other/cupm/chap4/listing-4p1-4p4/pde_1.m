  function ut=pde_1(t,u)
%
% Function pde_1 defines the ODEs in the MOL solution of the
% four-section retinal O2 transport model equations
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
% One vector to four vectors
  for i=1:nir uir(i)=u(i);             end
  for i=1:nor uor(i)=u(i+nir);         end
  for i=1:nfl ufl(i)=u(i+nir+nor);     end
  for i=1:ncc ucc(i)=u(i+nir+nor+nfl); end
%
% First order spatial derivatives
  uir_z=dss004(0,zl_ir,nir,uir); 
  uor_z=dss004(0,zl_or,nor,uor); 
  ufl_z=dss004(0,zl_fl,nfl,ufl); 
  ucc_z=dss004(0,zl_cc,ncc,ucc);   
%
% BCs
%
%   Inner retina
    uir(1)=pir_s;     
    uir(nir)=kir_or*uor(1);
%
%   Outer retina
    uor_z(1)=(Dir/Dor)*uir_z(nir);
    uor(nor)=kor_fl*ufl(1); 
%
%   Fluid layer
    ufl_z(1)=(Dor/Dfl)*uor_z(nor);
    ufl(nfl)=kfl_cc*ucc(1); 
%
%   Choriocapillaris
    ucc_z(1)=(Dfl/Dcc)*ufl_z(nfl);
    ucc(ncc)=pcc_s;      
%
% Second order spatial derivatives
  nl=1; nu=1;
  uir_zz=dss044(0,zl_ir,nir,uir,uir_z,nl,nu);
  nl=2; nu=1;  
  uor_zz=dss044(0,zl_or,nor,uor,uor_z,nl,nu);
  ufl_zz=dss044(0,zl_fl,nfl,ufl,ufl_z,nl,nu);
  ucc_zz=dss044(0,zl_cc,ncc,ucc,ucc_z,nl,nu);  
%
% PDEs
  uir_t=Dir*uir_zz-kir*uir;
  uor_t=Dor*uor_zz-kor*uor;
  ufl_t=Dfl*ufl_zz-kfl*ufl;  
  ucc_t=Dcc*ucc_zz-kcc*ucc; 
  uir_t(1)  =0;
  ucc_t(ncc)=0;  
%
% Four vectors to one vector
  for i=1:nir ut(i)            =uir_t(i); end
  for i=1:nor ut(i+nir)        =uor_t(i); end
  for i=1:nfl ut(i+nir+nor)    =ufl_t(i); end  
  for i=1:ncc ut(i+nir+nor+nfl)=ucc_t(i); end   
  ut=ut';
%
% Increment calls to pde_1
  ncall=ncall+1;
