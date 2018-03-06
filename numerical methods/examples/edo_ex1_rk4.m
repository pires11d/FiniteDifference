%Dados:
t0 = 0; h = 0.1; tf = 1; y0 = -1;

%Método de Runge-Kutta de 4a ordem:
t = t0;
y = y0;
tab = [t y];
while (t<tf)
	k1 = h*f(t,y); 
	k2 = h*f(t+h/2,y+k1/2);
	k3 = h*f(t+h/2,y+k2/2);
	k4 = h*f(t+h,y+k3);
	y = y + (k1+2*k2+2*k3+k4)/6;
	t = t + h;
        tab = [tab; t y];	
end
tab
tt = tab(:,1);
yy = tab(:,2);
y1 = -3*exp(-tt) - 2*(tt) + 2;
y1
plot(tt,yy,'o',tt,y1);

