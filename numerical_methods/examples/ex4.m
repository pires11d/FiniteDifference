%Dados
x = 1:0.1:3;
y = 1./(sqrt((4-x.^2).^2+0.02));
erro = 0.01*randn(size(x));
y = y.*(1+erro);

%Dados transformados
X = x.^2;
Y = 1./y.^2;
plot(X,Y,'o'), pause

%Ajuste
c = polyfit(x,y,2);

%Plot v2
Y = polyval(c,X);
	plot(X,Y), pause, close
