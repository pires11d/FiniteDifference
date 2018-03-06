clc
clf
clear all
close all

n = 20;
L = 0.001;
dx = L/(n-1);
x = 0:dx:L;

Da = 1e-6;
Db = 1e-6;
Dr = 1e-6;
k1 = 1;
k2 = 1;

Cao = 1.0;
Cbo = 0.9;
Cro = 0.5;

for i = 1:n
Ca(i) = 0;
Cb(i) = 0;
Cr(i) = 0;
Caa(i) = 0;
Cbb(i) = 0;
Crr(i) = 0;
end

Camin = 0;
Camax = 2;
Lmax = L;

Ca(1) = Cao;
Ca(2) = Ca(1);
Cb(n) = Cbo;
Cb(n-1) = Cb(n);
Cr(n) = Cro;
Cr(n-1) = Cr(n);

k = 1;
erro = 1;

while k<10
  for i=2:n-2
    Ca(i+1) = -(dx^2/Da)*(k1*Ca(i)*Cbb(i)+k2*Ca(i)*Crr(i))+2*Ca(i)-Ca(i-1);
    Ca(n)=Ca(n-1);
  end
  Caa = Ca;
  for i=n:-1:4
    Cb(i-2) = -(dx^2/Db)*(k1*Caa(i)*Cb(i))+2*Cb(i-1)-Cb(i);
    Cb(1)=Cb(2);
  end
  Cbb = Cb;
  for i=n:-1:4
    Cr(i-2) = -(dx^2/Dr)*(k2*Caa(i)*Cr(i)-k1*Caa(i)*Cbb(i))+2*Cr(i-1)-Cr(i);
    Cr(1)=Cr(2);
  end
  Crr = Cr;
  k = k+1;
end  


figure(1)
title('C = C(x)','FontSize',14)
plot(x,Caa,x,Cbb,x,Crr)
xlabel('x')
ylabel('C')
legend('Ca','Cb','Cr')
xlim([0 Lmax])
ylim([Camin,Camax])



