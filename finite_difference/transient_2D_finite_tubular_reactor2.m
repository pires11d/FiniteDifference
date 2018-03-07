clc
clf
clear all
close all

%% Malha
m = 9;          % Numero de pontos p/ z
n = 6;          % Numero de pontos p/ r

L = 3;          % Comprimento final
R = 0.5;        % Raio final
tf = 5;         % Tempo final

dz = L/(m-1);   
dr = R/(n-1);
dt = 0.1;

z = 0:dz:L;    
r = 0:dr:R;
t = 0:dt:tf;

mf = length(z);
nf = length(r);
pf = length(t);

%% aA + bB -> cC + dD
a = 1;
b = 1/2;
c = 1;
d = 1/2;
%% Constantes Cinéticas
ko = 20;
Ea = 10000;
R = 8.314;
dHr = -200000;
%% Constantes do Projeto
eps = 0.8;
u = 1.5;
dp = 0.01;
To = 500;
Te = To;
Tc = To;
Po = 1e5;
Pe = 2e5;
Dz = 0.0005;
Dr = 0.001;
rho = 1;
mu = 1e-3;
Cp = 4187;
cond = 1;
U = 10;
lambda_z = Dz*cond;
lambda_r = Dr*cond;
alpha_z = lambda_z/(rho*Cp);
alpha_r = lambda_r/(rho*Cp);
Re = rho*u*dp/mu;
f = 6.8*(((1-eps)^1.2)/(eps^3))*Re^(-0.2);

%% A
Ka = 1;
Cao = 0;
Cae = 1;
%% B 
Kb = 1;
Cbo = 0;
Cbe = 0.5;
%% C
Kc = 1;
Cco = 0;
Cce = 0;
%% D
Kd = 1;
Cdo = 0;
Cde = 0;

%% Limites
Ca = zeros(mf,nf,pf);
Cb = zeros(mf,nf,pf);
Cc = zeros(mf,nf,pf);
Cd = zeros(mf,nf,pf);
T = zeros(mf,nf,pf);
P = zeros(mf,pf);

Camin = 0;
Camax = 1;
Tmin = To;
Tmax = To+20;
Pmin = 0;
Pmax = Pe;

%% Condições Iniciais
Ca(:,:,1) = Cao;        %Ca(z,r,0) 
Cb(:,:,1) = Cbo;        %Cb(z,r,0)
Cc(:,:,1) = Cco;        %Cc(z,r,0)
Cd(:,:,1) = Cdo;        %Cd(z,r,0)
T(:,:,1) = To;          %T(z,r,0)
P(:,1) = Po;            %P(z,0)
%% Condições de Contorno na Borda Inicial (z = 0)
Ca(1,:,:) = Cae;        %Ca(0,r,t)
Cb(1,:,:) = Cbe;        %Cb(0,r,t)
Cc(1,:,:) = Cce;        %Cc(0,r,t)
Cd(1,:,:) = Cde;        %Cd(0,r,t)
T(1,:,:) = Te;          %T(0,r,t)
P(1,:) = Pe;            %P(0,t)            

%% DIFERENÇAS FINITAS    
for p = 1:pf-1
    for m = 2:mf-1
        for n = 2:nf-1
            k = ko*exp(-Ea/(R*T(m,n)));
            %% Condições de Contorno no Centro (r = 0)
            Ca(m,1,p) = 2*Ca(m,2,p)-Ca(m,3,p);
            Cb(m,1,p) = 2*Cb(m,2,p)-Cb(m,3,p);
            Cc(m,1,p) = 2*Cc(m,2,p)-Cc(m,3,p);
            Cd(m,1,p) = 2*Cd(m,2,p)-Cd(m,3,p);
            T(m,1,p) = 2*T(m,2,p)-T(m,3,p);     
            %% Balanço da espécie A
            adveccao_z_a = -u*(Ca(m,n,p)-Ca(m-1,n,p))/(dz);
            difusao_z_a = Dz*(Ca(m+1,n,p)-2*Ca(m,n,p)+Ca(m-1,n,p))/(dz^2);
            difusao_r_a = Dr*((1/((n-1)*(dr^2)))*(Ca(m,n+1,p)-Ca(m,n,p))+(Ca(m,n+1,p)-2*Ca(m,n,p)+Ca(m,n-1,p))/(dr^2));
            reacao_a = -(k/eps)*Ca(m,n,p)/(1+Ka*Ca(m,n,p)+Kb*Cb(m,n,p)+Kc*(Cc(m,n,p)+Kd*Cd(m,n,p)));
            Ca(m,n,p+1) = Ca(m,n,p) + dt*(adveccao_z_a + difusao_z_a + difusao_r_a + reacao_a);
            %% Balanço da espécie B:
            adveccao_z_b = -u*(Cb(m,n,p)-Cb(m-1,n,p))/(dz);
            difusao_z_b = Dz*(Cb(m+1,n,p)-2*Cb(m,n,p)+Cb(m-1,n,p))/(dz^2);
            difusao_r_b = Dr*((1/((n-1)*(dr^2)))*(Cb(m,n+1,p)-Cb(m,n,p))+(Cb(m,n+1,p)-2*Cb(m,n,p)+Cb(m,n-1,p))/(dr^2));
            reacao_b = reacao_a*(b/a);
            Cb(m,n,p+1) = Cb(m,n,p) + dt*(adveccao_z_b + difusao_z_b + difusao_r_b + reacao_b);
            %% Balanço da espécie C
            adveccao_z_c = -u*(Cc(m,n,p)-Cc(m-1,n,p))/(dz);
            difusao_z_c = Dz*(Cc(m+1,n,p)-2*Cc(m,n,p)+Cc(m-1,n,p))/(dz^2);
            difusao_r_c = Dr*((1/((n-1)*(dr^2)))*(Cc(m,n+1,p)-Cc(m,n,p))+(Cc(m,n+1,p)-2*Cc(m,n,p)+Cc(m,n-1,p))/(dr^2));
            reacao_c = -reacao_a*(c/a);
            Cc(m,n,p+1) = Cc(m,n,p) + dt*(adveccao_z_c + difusao_z_c + difusao_r_c + reacao_c);
            %% Balanço da espécie D
            adveccao_z_d = -u*(Cd(m,n,p)-Cd(m-1,n,p))/(dz);
            difusao_z_d = Dz*(Cd(m+1,n,p)-2*Cd(m,n,p)+Cd(m-1,n,p))/(dz^2);
            difusao_r_d = Dr*((1/((n-1)*(dr^2)))*(Cd(m,n+1,p)-Cd(m,n,p))+(Cd(m,n+1,p)-2*Cd(m,n,p)+Cd(m,n-1,p))/(dr^2));
            reacao_d = -reacao_a*(d/a);
            Cd(m,n,p+1) = Cd(m,n,p) + dt*(adveccao_z_d + difusao_z_d + difusao_r_d + reacao_d);     
            %% Balanço de Energia
            adveccao_z_T = -u*(T(m,n,p)-T(m-1,n,p))/(dz);
            conducao_z_T = alpha_z*(T(m+1,n,p)-2*T(m,n,p)+T(m-1,n,p))/(dz^2);
            conducao_r_T = alpha_r*((1/((n-1)*(dr^2)))*(T(m,n+1,p)-T(m,n,p))+(T(m,n+1,p)-2*T(m,n,p)+T(m,n-1,p))/(dr^2));
            reacao_T = reacao_a*dHr/(rho*Cp);
            troca_T = 2*(U/R)*(Tc - T(m,n,p));
            T(m,n,p+1) = T(m,n,p) + dt*(adveccao_z_T + conducao_z_T + conducao_r_T + reacao_T + troca_T);     
            %% Perda de Carga;
            adveccao_z_P = -u*(P(m,p)-P(m-1,p));
            atrito_P = -(f*rho*u^3)/dp;
            P(m,p+1) = P(m,p) + dt*(adveccao_z_P + atrito_P); 
        end
    end
    %% Condições de Contorno na Borda Final (z = L)
    Ca(mf,:,p)=2*Ca(mf-1,:,p)-Ca(mf-2,:,p);
    Cb(mf,:,p)=2*Cb(mf-1,:,p)-Cb(mf-2,:,p);
    Cc(mf,:,p)=2*Cc(mf-1,:,p)-Cc(mf-2,:,p);
    Cd(mf,:,p)=2*Cd(mf-1,:,p)-Cd(mf-2,:,p);
    T(mf,:,p)=2*T(mf-1,:,p)-T(mf-2,:,p);
    P(mf,p)=2*P(mf-1,p)-P(mf-2,p);
    %% Condições de Contorno na Parede (r = R)
    Ca(:,nf,p)=2*Ca(:,nf-1,p)-Ca(:,nf-2,p);
    Cb(:,nf,p)=2*Cb(:,nf-1,p)-Cb(:,nf-2,p);
    Cc(:,nf,p)=2*Cc(:,nf-1,p)-Cc(:,nf-2,p);
    Cd(:,nf,p)=2*Cd(:,nf-1,p)-Cd(:,nf-2,p);
    T(:,nf,p)=2*T(:,nf-1,p)-T(:,nf-2,p);
end

%% Gráfico Dinâmico das Concentrações Ci(z,r,t)
for p = 1:pf-1    
figure(1)

subplot(2,2,1)
surf(r,z,Ca(:,:,p))
hold on
title('Concentração de A','FontSize',14)
surf(-r,z,Ca(:,:,p))
xlabel('R raio (m)')
ylabel('Z comprimento (m)')
zlabel('Ca (mol/L)')
caxis([Camin,Camax])
view(-100,20)
hold off
xlim([-L L])
ylim([0 L])
zlim([Camin Camax])

subplot(2,2,2)
surf(r,z,Cb(:,:,p))
hold on
surf(-r,z,Cb(:,:,p))
title('Concentração de B','FontSize',14)
xlabel('R raio (m)')
ylabel('Z comprimento (m)')
zlabel('Cb (mol/L)')
caxis([Camin,Camax])
view(-100,20)
hold off
xlim([-L L])
ylim([0 L])
zlim([Camin Camax])

subplot(2,2,3)
surf(r,z,Cc(:,:,p))
hold on
surf(-r,z,Cc(:,:,p))
title('Concentração de C','FontSize',14)
xlabel('R raio (m)')
ylabel('Z comprimento (m)')
zlabel('Cc (mol/L)')
caxis([Camin,Camax])
view(-80,20)
hold off
xlim([-L L])
ylim([0 L])
zlim([Camin Camax])

subplot(2,2,4)
surf(r,z,Cd(:,:,p))
hold on
surf(-r,z,Cd(:,:,p))
title('Concentração de D','FontSize',14)
xlabel('R raio (m)')
ylabel('Z comprimento (m)')
zlabel('Cd (mol/L)')
caxis([Camin,Camax])
view(-80,20)
hold off
xlim([-L L])
ylim([0 L])
zlim([Camin Camax])
drawnow
end

%% Gráfico Dinâmico da Temperatura T(z,r,t)
for p=1:pf-1
figure(2)

subplot(1,1,1)
surf(r,z,T(:,:,p))
hold on
title('Temperatura ao longo do reator','FontSize',14)
surf(-r,z,T(:,:,p))
xlabel('R raio (m)')
ylabel('Z comprimento (m)')
zlabel('T (K)')
caxis([Tmin,Tmax])
view(-100,20)
hold off
xlim([-L L])
ylim([0 L])
zlim([Tmin Tmax])
drawnow
end

%% Gráfico Dinâmico da Pressão P(z,r,t)
for p=1:pf-1
figure(3)

subplot(1,1,1)
plot(z,P(:,p)/100000)
title('Pressão ao longo do reator','FontSize',14)
xlabel('Z comprimento (m)')
ylabel('P (bar)')
xlim([0 L])
ylim([Pmin Pmax]/100000)
drawnow
end


