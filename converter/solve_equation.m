pkg load symbolic
clc
more("off")

% Definição de Variáveis
syms uc u A d rho mar mox m

% Resolve a Equação Simbólica
A = pi()*d^2/4
uc = 1.4*(sqrt(d)*(mox/mar)^0.25)
u = uc+5
m = rho*A*u

% Resolve a Equação com Valores
solve(m-mar-mox==0,mar)

