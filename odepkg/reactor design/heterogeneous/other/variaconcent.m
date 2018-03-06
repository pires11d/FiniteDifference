%Samara Silva de Souza
%email: samara@intelab.ufsc.br
%Exercício 03/09/2014
%Conjunto de equações diferenciais que descreve a mudança da concentração
%das espécies A, B e C em função do tempo.

function variaconcent

clc;
clear;
global k1 k2

% tspan representa o intervalo de integração no tempo (ou var.
% independente)
tspan=[0 5];
% c0 representa o vetor de valores iniciais da var. dependente, nesse caso
% concentração.
c0=[5 0 0];

[t,c]=ode45(@concent,tspan,c0);

plot(t,c(:,1),'+g',t,c(:,2),'b*',t,c(:,3),'r')
title('Concentração vs tempo');
xlabel('t(h)');
ylabel('Concentração de cada espécie');  
legend('Ca','Cb','Cc');

function dcdt= concent(t,c)
% .m file que define a função a ser integrada. dcdt representa a equação a
% ser integrada, e c é a variável dependente do problema e t representa a
% relação da derivada representada na equação diferencial.

global k1 k2

k1=1;
k2=2;
dcdt=[-k1*c(1);k1*c(1)-k2*c(2); k2*c(2)];
