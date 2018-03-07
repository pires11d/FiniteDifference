%Dados
x = 1:0.05:3;
y = 1./(sqrt((4-x.^2).^2+0.02));
erro = 0.01*randn(size(x));
y = y.*(1+erro);
plot(x,y,'o'), pause, hold on

%Ajuste
c4 = polyfit(x,y,4);
c8 = polyfit(x,y,8);
c12 = polyfit(x,y,12);

%Plot
xx = 1:0.01:3;
yy4 = polyval(c4,xx);
	plot(xx,yy4), pause
yy8 = polyval(c8,xx);
	plot(xx,yy8), pause
yy12 = polyval(c12,xx);
	plot(xx,yy12), pause

%Dados transformados
X = x.^2;
Y = 1./y.^2;
plot(x,y,'o'), pause

%Ajuste v2
c = polyfit(x,y,2);

%Plot v2
yy = polyval(c,xx);
	plot(xx,yy), pause, close
