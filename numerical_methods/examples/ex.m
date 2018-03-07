function f = ex (t,x)
k1=1;
k2=0;
k3=2;
k4=3;
f = zeros(3,1);
f(1)=-k1*x(1)+k2*x(2);
f(2)=k1*x(1)-k2*x(2)-k3*x(2)+k4*x(3);
f(3)=k3*x(2)-k4*x(3);
end
