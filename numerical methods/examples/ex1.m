%Dados
n = 21;
x = linspace(0,4,n);
y = exp(x);
plot (x,y,'o'), hold on


%Ajuste
A = [	n sum(x) sum(x.^2) sum(x.^3); 
	sum(x)	sum(x.^2) sum(x.^3) sum(x.^4); 
	sum(x.^2) sum(x.^3) sum(x.^4) sum(x.^5)
	sum(x.^3) sum(x.^4) sum(x.^5) sum(x.^6)];
b = [	sum(y); 
	sum(x.*y); 
	sum(x.^2.*y)
	sum(x.^3.*y)];
c = A\b
cond(A)

%Plot
yy = c(1) + c(2)*x + c(3)*x.^2 + c(4)*x.^3;
plot (x,yy)
