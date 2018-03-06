  function S=jpattern_num
%
  global n
%
% Sparsity pattern of the Jacobian matrix based on a
% numerical evaluation.  Note: the reference to the ODE
% routine (two places below) should be edited to specify
% the current ODE routine.
%
% Set independent, dependent variables for the calculation 
% of the sparsity pattern
  tbase=0;
  for i=1:n
    ybase(i)=0.5;
  end
  ybase=ybase';
%
% Compute the corresponding derivative vector
  ytbase=pde_1(tbase,ybase);
  fac=[];
  thresh=1e-16;
  vectorized='on';
  [Jac,fac]=numjac(@pde_1,tbase,ybase,ytbase,thresh,fac,vectorized);
%
% Replace nonzero elements by "1" 
  S=spones(sparse(Jac));
%
% Plot the map
  figure(3)
  spy(S);
  xlabel('dependent variables');
  ylabel('semi-discrete equations');
%
% Compute the percentage of non-zero elements
  [njac,mjac]=size(S);
  ntotjac=njac*mjac;
  non_zero=nnz(S);
  non_zero_percent=non_zero/ntotjac*100;
  stat=sprintf('Jacobian sparsity pattern - nonzeros %d (%.3f%%)',...
       non_zero,non_zero_percent);
  title(stat);