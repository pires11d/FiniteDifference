%Samara Silva de Souza
%email: samara@intelab.ufsc.br
%Exerc�cio 03/09/2014
%Conjunto de equa��es diferenciais que descreve a mudan�a da concentra��o
%das esp�cies A, B e C em fun��o do tempo.

function variaconcent

clc;
clear;
global k1 k2

% tspan representa o intervalo de integra��o no tempo (ou var.
% independente)
tspan=[0 5];
% c0 representa o vetor de valores iniciais da var. dependente, nesse caso
% concentra��o.
c0=[5 0 0];

[t,c]=ode45(@concent,tspan,c0);

plot(t,c(:,1),'+g',t,c(:,2),'b*',t,c(:,3),'r')
title('Concentra��o vs tempo');
xlabel('t(h)');
ylabel('Concentra��o de cada esp�cie');  
legend('Ca','Cb','Cc');

function dcdt= concent(t,c)
% .m file que define a fun��o a ser integrada. dcdt representa a equa��o a
% ser integrada, e c � a vari�vel dependente do problema e t representa a
% rela��o da derivada representada na equa��o diferencial.

global k1 k2

k1=1;
k2=2;
dcdt=[-k1*c(1);k1*c(1)-k2*c(2); k2*c(2)];
