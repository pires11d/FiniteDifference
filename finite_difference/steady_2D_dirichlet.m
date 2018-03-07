clf
clc
close all
clear all

tol = 10e-5;
kmax = 20;
m = 3;

T0j = 75;
Ti0 = 0;
T4j = 50;
Ti4 = 100;

A= [-4 1 0 1 0 0 0 0 0
    1 -4 1 0 1 0 0 0 0
    0 1 -4 0 0 1 0 0 0
    1 0 0 -4 1 0 1 0 0
    0 1 0 1 -4 1 0 1 0
    0 0 1 0 1 -4 0 0 1
    0 0 0 1 0 0 -4 1 0
    0 0 0 0 1 0 1 -4 1
    0 0 0 0 0 1 0 1 -4];
b= [-75
    -75
    -175
    0
    0
    -100
    -50
    -50
    -150];

TT = EG(A,b);

% for i=2:n+1
%     for j=2:n+1
%         T(i,j)=TT 
%         
%         for i=2:n+1

T = [(Ti4+T0j)/2    Ti4     Ti4     Ti4        (Ti4+T4j)/2
    T0j   TT(3,1) TT(6,1) TT(9,1)  T4j
    T0j   TT(2,1) TT(5,1) TT(8,1)  T4j
    T0j   TT(1,1) TT(4,1) TT(7,1)  T4j
    (T0j+Ti0)/2     Ti0     Ti0     Ti0        (Ti0+T4j)/2]

x = linspace(0,1,5);
y = linspace(0,1,5);
xx = linspace(0,1,3);
yy = linspace(0,1,3);

figure(1)
surf(x,y,T)
xlabel('x')
ylabel('y')
zlabel('T')

k = 0.49;
dx = 10;
dy = 10;

qx = zeros(m,m);
n = sqrt(length(A));

for j=2:n+1
    for i=2:n+1
        qx(i-1,j-1) = (k/(2*dx))*(T(i+1,j)-T(i-1,j));
    end
end

figure(2)
qx
surf(xx,yy,qx)
xlabel('x')
ylabel('y')
zlabel('qx')

for j=2:n+1
    for i=2:n+1
        qy(i-1,j-1) = (k/(2*dy))*(T(i,j+1)-T(i,j-1));
    end
end

figure(3)
qy
surf(xx,yy,qy)
xlabel('x')
ylabel('y')
zlabel('qy')

figure(4)
qn = sqrt(qx.^2+qy.^2)
surf(xx,yy,qn)
xlabel('x')
ylabel('y')
zlabel('qn')


if qx >= 0
   theta = atan(qy./qx)*180/pi;
else
   theta = atan(qy./qx)*pi*180/pi;
end
