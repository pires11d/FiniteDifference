clf
clc
close all
clear all

tol = 10e-5;
kmax = 20;

T00 = 0;
T0j = 75;
Ti0 = 0;
T4j = 50;
Ti4 = 100;

k = 0.49;
dx = 10;
dy = 10;

A= [-4 2 0 0 1 0 0 0 0 0 0 0
    1 -4 1 0 0 1 0 0 0 0 0 0
    0 1 -4 1 0 0 1 0 0 0 0 0
    0 0 1 -4 0 0 0 1 0 0 0 0 
    1 0 0 0 -4 2 0 0 1 0 0 0
    0 1 0 0 1 -4 1 0 0 1 0 0
    0 0 1 0 0 1 -4 1 0 0 1 0
    0 0 0 1 0 0 1 -4 0 0 0 1
    0 0 0 0 1 0 0 0 -4 2 0 0
    0 0 0 0 0 1 0 0 1 -4 1 0
    0 0 0 0 0 0 1 0 0 1 -4 1
    0 0 0 0 0 0 0 1 0 0 1 -4]

b= [-T00
    -T0j
    -T0j
    -T0j-Ti4
    0
    0
    0
    -Ti4
    -T4j
    -T4j
    -T4j
    -T4j-Ti4];

TT = EG(A,b);

% for i=2:n+1
%     for j=2:n+1
%         T(i,j)=TT 
%         
%         for i=2:n+1

T = [(Ti4+T0j)/2  Ti4     Ti4     Ti4       (Ti4+T4j)/2
    T0j   TT(4,1) TT(8,1) TT(12,1)  T4j
    T0j   TT(3,1) TT(7,1) TT(11,1)  T4j
    T0j   TT(2,1) TT(6,1) TT(10,1)  T4j
    (T0j+TT(1,1))/2   TT(1,1) TT(5,1) TT(9,1)   (T4j+TT(9,1))/2]

x = linspace(0,1,5);
y = linspace(0,1,5);
surf(x,y,T)
xlabel('x')
ylabel('y')
zlabel('T')

qx = zeros(3,3);
n = sqrt(length(A));

for i=2:n+1
    for j=2:n+1
        qx(i-1,j-1) = -(k/(2*dx))*(T(i+1,j)-T(i-1,j));
    end
end

qx

for j=2:n+1
    for i=2:n+1
        qy(i-1,j-1) = -(k/(2*dy))*(T(i,j+1)-T(i,j-1));
    end
end

qy

qn = sqrt(qx.^2+qy.^2)

if qx >= 0
   theta = atan(qy./qx)*180/pi
else
   theta = atan(qy./qx)*pi*180/pi
end
