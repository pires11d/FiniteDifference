function [raiz,k] = secante (x0,nrep,tol)

fx0 = f(x0);
x1 = x0 + 1*e-3;
erro = 2*tol;
k = 0;

while (k < nrep && erro > tol)
fx1 = f(x1);
x2 = x1 - fx1*((x1-x0)/(fx1-fx0));
erro = abs(x2-x1);
x1 = x2;
k = k+1;
end

raiz = x1;
end
