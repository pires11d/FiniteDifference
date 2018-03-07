%Dados
x = 0:0.5:6;
y = -x.^3 + 6*x.^2 + 2;
erro = randn(size(x));
y = y+erro;
plot(x,y,'o'), hold on

%Ajuste
c = polyfit(x,y,12)

%Plot
xx = 0:0.01:6;
yy = polyval(c,xx);
plot(xx,yy)
