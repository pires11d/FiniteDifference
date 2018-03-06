clear all
global k1 k2
k1=1; k2=2;

% tspan representa o intervalo de integração no tempo (ou var.
% independente)
tspan=[0 5];
% c0 representa o vetor de valores iniciais da var. dependente, nesse caso
% concentração.
c0=[5 0 0];

[t,c]=ode45(@Ex32,tspan,c0);

plot(t,c(:,1),'+g',t,c(:,2),'b*',t,c(:,3),'r')
title('Concentração vs tempo');
xlabel('t(h)');
ylabel('Concentração de cada espécie');  
legend('Ca','Cb','Cc');

