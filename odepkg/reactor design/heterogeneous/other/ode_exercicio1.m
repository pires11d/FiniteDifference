clear all
global k1 k2
k1=1; k2=2;

% tspan representa o intervalo de integra��o no tempo (ou var.
% independente)
tspan=[0 5];
% c0 representa o vetor de valores iniciais da var. dependente, nesse caso
% concentra��o.
c0=[5 0 0];

[t,c]=ode45(@Ex32,tspan,c0);

plot(t,c(:,1),'+g',t,c(:,2),'b*',t,c(:,3),'r')
title('Concentra��o vs tempo');
xlabel('t(h)');
ylabel('Concentra��o de cada esp�cie');  
legend('Ca','Cb','Cc');

