%Dados:
t0 = 0; h = 0.1; tf = 1; y0 = -1;

%Método de Euler:
t = t0;
y = y0;
tab = [t y];
while (t<tf)
	y = y + h*f(t,y);
	t = t + h;
	tab = [tab; t y];
end
tab
tt = tab(:,1);
yy = tab(:,2);
y1 = -3*e.^-(tt) - 2*(tt) + 2;
plot (tt,yy,tt,y1)

