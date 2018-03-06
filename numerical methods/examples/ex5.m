%Dados
x = [-1 0 1];
y = [3.1 0.9 2.9];
plot (x,y,'o'), hold on


%Ajuste
A = [	3 sum(x) sum(x.^2); 
	sum(x)	sum(x.^2) sum(x.^3); 
	sum(x.^2) sum(x.^3) sum(x.^4)];
b = [	sum(y); 
	sum(x.*y); 
	sum(x.^2.*y)];
c = A\b

%Plot
xx = [-10:0.01:10];
yy = c(1)*xx.^2 + c(2)*xx + c(3);
plot (xx,yy)


