pkg load symbolic
clc
more('off')
warning('off', 'all');

% Abre arquivo com Conversão de Unidades
unit_conversion

% Definição de Variáveis
syms u A

% Conversão de Unidades
u = u*_ft/_min;
A = A*_ft2;

% Resolve a Equação Simbólica
Q = u*A;
vpa(Q,4)

warning('on', 'all');
