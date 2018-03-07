clc
clf
clear all
close all

%% Malha
m = 8;          % Numero de pontos p/ z
n = 4;          % Numero de pontos p/ r
q = 4;          % Numero de pontos p/ rp

L = 3;          % Comprimento final
R = 0.5;        % Raio final
tf = 5;         % Tempo final
Rp = R/10;      % Raio final da partícula

dz = L/(m-1);   
dr = R/(n-1);
dt = 0.1;
drp = Rp/(q-1);

z = 0:dz:L;    
r = 0:dr:R;
t = 0:dt:tf;
rp = 0:drp:Rp;

mf = length(z);
nf = length(r);
pf = length(t);
qf = length(rp);

%% aA + bB -> cC + dD
a = 1;
b = 1/2;
c = 1;
d = 1/2;
%% Constantes Cinéticas
ko = 10;
Ea = 10000;
R = 8.314;
dHr = -20000;
%% Constantes do Projeto
eps = 0.4;
u = 1.2;
dp = 2*Rp;
eps_p = 0.8;
To = 300;
Te = To;
Tc = To;
dT = 20;
Po = 1e5;
Pe = 2e5;
Dz = 0.001;
Dr = 0.001;
Dr_p = 0.0001;
rho = 5;
rho_p = 5000;
mu = 1e-5;
Cp = 10;
Cp_p = 1;
cond = 10;
cond_p = 100;
U = 20;
lambda_z = Dz*cond;
lambda_r = Dr*cond;
lambda_rp = Dr_p*cond_p;
alpha_z = lambda_z/(rho*Cp);
alpha_r = lambda_r/(rho*Cp);
alpha_rp = lambda_rp/(rho_p*Cp_p);
g = 9.81;
Re = rho*u*dp/mu;
f = 6.8*(((1-eps)^1.2)/(eps^3))*Re^(-0.2);
eta = 1;

%% A
Ka = 0;
Cao = 0;
Cae = 1;
%% B 
Kb = 0;
Cbo = 0;
Cbe = 0.5;
%% C
Kc = 0;
Cco = 0;
Cce = 0;
%% D
Kd = 0;
Cdo = 0;
Cde = 0;

%% Limites
Ca = zeros(mf,nf,pf);
Cb = zeros(mf,nf,pf);
Cc = zeros(mf,nf,pf);
Cd = zeros(mf,nf,pf);
T = zeros(mf,nf,pf);
P = zeros(mf,pf);
Cap = zeros(mf,qf,pf);
Cbp = zeros(mf,qf,pf);
Ccp = zeros(mf,qf,pf);
Cdp = zeros(mf,qf,pf);
Tp = zeros(mf,qf,pf);

Camin = 0;
Camax = 1;
Tmin = To-2*dT;
Tmax = To+2*dT;
Tmin2 = To-dT;
Tmax2 = To+dT;
Pmin = 0;
Pmax = Pe;

%% Condições Iniciais
Ca(:,:,1) = Cao;        %Ca(z,r,0) 
Cb(:,:,1) = Cbo;        %Cb(z,r,0)
Cc(:,:,1) = Cco;        %Cc(z,r,0)
Cd(:,:,1) = Cdo;        %Cd(z,r,0)
T(:,:,1) = To;          %T(z,r,0)
P(:,1) = Po;            %P(z,0)
Cap(:,:,1) = Cao;       %Cap(z,rp,0)
Cbp(:,:,1) = Cao;       %Cbp(z,rp,0)
Ccp(:,:,1) = Cao;       %Ccp(z,rp,0)
Cdp(:,:,1) = Cao;       %Cdp(z,rp,0)
Tp(:,:,1) = To;         %Tp(z,rp,0)

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
            for q = 2:qf-1
            k = ko*exp(-Ea/(R*T(m,n,p)));
            kp = ko*exp(-Ea/(R*Tp(n,q,p)));
            %% Condições de Contorno no Centro (r = 0)
            Ca(m,1,p)=2*Ca(m,2,p)-Ca(m,3,p);
            Cb(m,1,p)=2*Cb(m,2,p)-Cb(m,3,p);
            Cc(m,1,p)=2*Cc(m,2,p)-Cc(m,3,p);
            Cd(m,1,p)=2*Cd(m,2,p)-Cd(m,3,p);
            T(m,1,p)=2*T(m,2,p)-T(m,3,p);     
            %% Condições de Contorno no Pellet
            Cap(m,1,p) = 2*Cap(m,2,p)-Cap(m,3,p);
            Cap(m,qf,p) = Ca(m,1,p);
            Cap(1,:,p) = Ca(1,1,p);
            Cbp(m,1,p) = 2*Cbp(m,2,p)-Cbp(m,3,p);
            Cbp(m,qf,p) = Cb(m,1,p);
            Cbp(1,:,p) = Cb(1,1,p);
            Ccp(m,1,p) = 2*Ccp(m,2,p)-Ccp(m,3,p);
            Ccp(m,qf,p) = Cc(m,1,p);
            Ccp(1,:,p) = Cc(1,1,p);
            Cdp(m,1,p) = 2*Cdp(m,2,p)-Cdp(m,3,p);
            Cdp(m,qf,p) = Cd(m,1,p);
            Cdp(1,:,p) = Cd(1,1,p);
            Tp(m,1,p) = 2*Tp(m,2,p)-Tp(m,3,p);
            Tp(m,qf,p) = T(m,1,p);
            Tp(1,:,p) = T(1,1,p);
            %% Balanço da espécie A:
            adveccao_z_a = -u*(Ca(m,n,p)-Ca(m-1,n,p))/(dz);
            difusao_z_a = Dz*(Ca(m+1,n,p)-2*Ca(m,n,p)+Ca(m-1,n,p))/(dz^2);
            difusao_r_a = Dr*((1/((n-1)*(dr^2)))*(Ca(m,n+1,p)-Ca(m,n,p))+(Ca(m,n+1,p)-2*Ca(m,n,p)+Ca(m,n-1,p))/(dr^2));
            reacao_a = -k*Ca(m,n,p)/(1+Ka*Ca(m,n,p)+Kb*Cb(m,n,p)+Kc*(Cc(m,n,p)+Kd*Cd(m,n,p)));
            Ca(m,n,p+1) = Ca(m,n,p) + dt*(adveccao_z_a + difusao_z_a + difusao_r_a + (eta/eps)*reacao_a);
            %% Balanço da espécie B
            adveccao_z_b = -u*(Cb(m,n,p)-Cb(m-1,n,p))/(dz);
            difusao_z_b = Dz*(Cb(m+1,n,p)-2*Cb(m,n,p)+Cb(m-1,n,p))/(dz^2);
            difusao_r_b = Dr*((1/((n-1)*(dr^2)))*(Cb(m,n+1,p)-Cb(m,n,p))+(Cb(m,n+1,p)-2*Cb(m,n,p)+Cb(m,n-1,p))/(dr^2));
            reacao_b = reacao_a*(b/a);
            Cb(m,n,p+1) = Cb(m,n,p) + dt*(adveccao_z_b + difusao_z_b + difusao_r_b + (eta/eps)*reacao_b);
            %% Balanço da espécie C
            adveccao_z_c = -u*(Cc(m,n,p)-Cc(m-1,n,p))/(dz);
            difusao_z_c = Dz*(Cc(m+1,n,p)-2*Cc(m,n,p)+Cc(m-1,n,p))/(dz^2);
            difusao_r_c = Dr*((1/((n-1)*(dr^2)))*(Cc(m,n+1,p)-Cc(m,n,p))+(Cc(m,n+1,p)-2*Cc(m,n,p)+Cc(m,n-1,p))/(dr^2));
            reacao_c = -reacao_a*(c/a);
            Cc(m,n,p+1) = Cc(m,n,p) + dt*(adveccao_z_c + difusao_z_c + difusao_r_c + (eta/eps)*reacao_c);
            %% Balanço da espécie D
            adveccao_z_d = -u*(Cd(m,n,p)-Cd(m-1,n,p))/(dz);
            difusao_z_d = Dz*(Cd(m+1,n,p)-2*Cd(m,n,p)+Cd(m-1,n,p))/(dz^2);
            difusao_r_d = Dr*((1/((n-1)*(dr^2)))*(Cd(m,n+1,p)-Cd(m,n,p))+(Cd(m,n+1,p)-2*Cd(m,n,p)+Cd(m,n-1,p))/(dr^2));
            reacao_d = -reacao_a*(d/a);
            Cd(m,n,p+1) = Cd(m,n,p) + dt*(adveccao_z_d + difusao_z_d + difusao_r_d + (eta/eps)*reacao_d);     
            %% Balanço de Energia
            adveccao_z_T = -u*(T(m,n,p)-T(m-1,n,p))/(dz);
            conducao_z_T = alpha_z*(T(m+1,n,p)-2*T(m,n,p)+T(m-1,n,p))/(dz^2);
            conducao_r_T = alpha_r*((1/((n-1)*(dr^2)))*(T(m,n+1,p)-T(m,n,p))+(T(m,n+1,p)-2*T(m,n,p)+T(m,n-1,p))/(dr^2));
            reacao_T = reacao_a*dHr/(rho*Cp);
            troca_T = 2*(U/R)*(Tc - T(m,n,p));
            T(m,n,p+1) = T(m,n,p) + dt*(adveccao_z_T + conducao_z_T + conducao_r_T + (eta/eps)*reacao_T + troca_T);     
            %% Perda de Carga
            adveccao_z_P = -u*(P(m,p)-P(m-1,p));
            atrito_P = -(f*rho*u^3)/(dp*g);
            P(m,p+1) = P(m,p) + dt*(adveccao_z_P + atrito_P); 
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Balanço de A na partícula
            difusao_rp_a = Dr_p*((Cap(m,q+1,p)-2*Cap(m,q,p)+Cap(m,q-1,p))/(drp^2)+(2/((q-1)*(drp^2)))*(Cap(m,q,p)-Cap(m,q-1,p)));
            reacao_p_a = -kp*Cap(m,q,p)/(1+Ka*Cap(m,q,p)+Kb*Cbp(m,q,p)+Kc*Ccp(m,q,p)+Kd*Cdp(m,q,p));
            Cap(m,q,p+1) = Cap(m,q,p) + dt*((1/eps_p)*(difusao_rp_a) + (1/eps_p)*reacao_p_a); 
            %% Balanço de B na partícula
            difusao_rp_b = Dr_p*((Cbp(m,q+1,p)-2*Cbp(m,q,p)+Cbp(m,q-1,p))/(drp^2)+(2/((q-1)*(drp^2)))*(Cbp(m,q,p)-Cbp(m,q-1,p)));
            reacao_p_b = reacao_p_a*(b/a);
            Cbp(m,q,p+1) = Cbp(m,q,p) + dt*((1/eps_p)*(difusao_rp_b) + (1/eps_p)*reacao_p_b); 
            %% Balanço de C na partícula
            difusao_rp_c = Dr_p*((Ccp(m,q+1,p)-2*Ccp(m,q,p)+Ccp(m,q-1,p))/(drp^2)+(2/((q-1)*(drp^2)))*(Ccp(m,q,p)-Ccp(m,q-1,p)));
            reacao_p_c = reacao_p_a*(-c/a);
            Ccp(m,q,p+1) = Ccp(m,q,p) + dt*((1/eps_p)*(difusao_rp_c) + (1/eps_p)*reacao_p_c); 
            %% Balanço de D na partícula
            difusao_rp_d = Dr_p*((Cdp(m,q+1,p)-2*Cdp(m,q,p)+Cdp(m,q-1,p))/(drp^2)+(2/((q-1)*(drp^2)))*(Cdp(m,q,p)-Cdp(m,q-1,p)));
            reacao_p_d = reacao_p_a*(-d/a);
            Cdp(m,q,p+1) = Cdp(m,q,p) + dt*((1/eps_p)*(difusao_rp_d) + (1/eps_p)*reacao_p_d); 
            %% Balanço de Energia na partícula
            conducao_rp_T = alpha_rp*((2/((q-1)*(drp^2)))*(Tp(m,q,p)-Tp(m,q-1,p))+(Tp(m,q+1,p)-2*Tp(m,q,p)+Tp(m,q-1,p))/(drp^2));
            reacao_Tp = reacao_p_a*dHr/(rho_p*Cp_p);
            Tp(m,q,p+1) = Tp(m,q,p) + dt*(conducao_rp_T + reacao_Tp);         
            end
        end
    end
    %% Condições de Contorno na Borda Final (z = L)
    Ca(mf,:,p)=2*Ca(mf-1,:,p)-Ca(mf-2,:,p);
    Cb(mf,:,p)=2*Cb(mf-1,:,p)-Cb(mf-2,:,p);
    Cc(mf,:,p)=2*Cc(mf-1,:,p)-Cc(mf-2,:,p);
    Cd(mf,:,p)=2*Cd(mf-1,:,p)-Cd(mf-2,:,p);
    T(mf,:,p)=2*T(mf-1,:,p)-T(mf-2,:,p);
    P(mf,p)=2*P(mf-1,p)-P(mf-2,p);
    Cap(mf,:,p)=2*Cap(mf-1,:,p)-Cap(mf-2,:,p);
    Cbp(mf,:,p)=2*Cbp(mf-1,:,p)-Cbp(mf-2,:,p);
    Ccp(mf,:,p)=2*Ccp(mf-1,:,p)-Ccp(mf-2,:,p);
    Cdp(mf,:,p)=2*Cdp(mf-1,:,p)-Cdp(mf-2,:,p);
    Tp(mf,:,p)=2*Tp(mf-1,:,p)-Tp(mf-2,:,p);
    %% Condições de Contorno na Parede (r = R)
    Ca(:,nf,p)=2*Ca(:,nf-1,p)-Ca(:,nf-2,p);
    Cb(:,nf,p)=2*Cb(:,nf-1,p)-Cb(:,nf-2,p);
    Cc(:,nf,p)=2*Cc(:,nf-1,p)-Cc(:,nf-2,p);
    Cd(:,nf,p)=2*Cd(:,nf-1,p)-Cd(:,nf-2,p);
    T(:,nf,p)=2*T(:,nf-1,p)-T(:,nf-2,p);
end

eta_local = ones(mf,qf);
reacao_p_a_SS = ones(mf,qf);
reacao_a_SS = ones(mf);

%% Cálculo do Eta
for q = 1:qf-1
   	for m = 2:mf-1
        reacao_p_a_SS(m,q) = -k*Cap(m,q,pf-1);
      	reacao_a_SS(m) = -k*Ca(m,2,pf-1);
      	eta_local(m,q) = reacao_p_a_SS(m,q)/reacao_a_SS(m);
  	end
end

%% Curva Dinâmica da Temperatura Tp (rp,t)
for p = 1:(pf-1)/2   
figure(1)

plot(rp,Tp(end,:,p))
hold on
title('Temperatura no Pellet em z = L','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Tp (K)')
xlim([0 Rp])
ylim([Tmin2 Tmax2])
hold off
drawnow
% set(gcf, 'Position', get(0,'Screensize'));
end

%% Curvas Dinâmicas das Concentrações Cip(rp,t)
for p = 1:(pf-1)/2    
figure(2)

subplot(2,2,1)
plot(rp,Cap(end,:,p))
hold on
title('Ca no Pellet em z = L','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Cap (mol/L)')
xlim([0 Rp])
ylim([Camin Camax])
hold off

subplot(2,2,2)
plot(rp,Cbp(end,:,p))
hold on
title('Cb no Pellet em z = L','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Cbp (mol/L)')
xlim([0 Rp])
ylim([Camin Camax])
hold off

subplot(2,2,3)
plot(rp,Ccp(end,:,p))
hold on
title('Cc no Pellet em z = L','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Ccp (mol/L)')
xlim([0 Rp])
ylim([Camin Camax])
hold off

subplot(2,2,4)
plot(rp,Cdp(end,:,p))
hold on
title('Cd no Pellet em z = L','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Cdp (mol/L)')
xlim([0 Rp])
ylim([Camin Camax])
hold off
drawnow
end

%% Superfícies Dinâmicas das Concentrações dentro do Pellet Cip(z,r,t)
for p = 1:(pf-1)/2   
figure(3)

subplot(2,2,1)
surf(rp,z(z~=0),Cap(z~=0,:,p))
hold on
title('Ca no pellet','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Z comprimento (m)')
zlabel('Cap (mol/L)')
caxis([Camin,Camax])
view(-100,20)
hold off
xlim([0 Rp])
ylim([0 L])
zlim([Camin Camax])

subplot(2,2,2)
surf(rp,z(z~=0),Cbp(z~=0,:,p))
hold on
title('Cb no pellet','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Z comprimento (m)')
zlabel('Cbp (mol/L)')
caxis([Camin,Camax])
view(-100,20)
hold off
xlim([0 Rp])
ylim([0 L])
zlim([Camin Camax])

subplot(2,2,3)
surf(rp,z(z~=0),Ccp(z~=0,:,p))
hold on
title('Cc no pellet','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Z comprimento (m)')
zlabel('Ccp (mol/L)')
caxis([Camin,Camax])
view(-80,20)
hold off
xlim([0 Rp])
ylim([0 L])
zlim([Camin Camax])

subplot(2,2,4)
surf(rp,z(z~=0),Cdp(z~=0,:,p))
hold on
title('Cd no pellet','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Z comprimento (m)')
zlabel('Cdp (mol/L)')
caxis([Camin,Camax])
view(-80,20)
hold off
xlim([0 Rp])
ylim([0 L])
zlim([Camin Camax])
drawnow
end

%% Superfícies Dinâmicas das Concentrações Ci(z,r,t)
for p = 1:(pf-1)/2    
figure(4)

subplot(2,2,1)
surf(r,z,Ca(:,:,p))
hold on
title('Concentração A no Reator','FontSize',14)
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
title('Concentração B no Reator','FontSize',14)
surf(-r,z,Cb(:,:,p))
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
title('Concentração C no Reator','FontSize',14)
surf(-r,z,Cc(:,:,p))
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
title('Concentração D no Reator','FontSize',14)
surf(-r,z,Cd(:,:,p))
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

%% Superfície Dinâmica da Temperatura T(z,r,t)
for p=1:(pf-1)/2
figure(5)

subplot(1,1,1)
surf(r,z,T(:,:,p))
hold on
title('Temperatura no Reator','FontSize',14)
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

%% Curva Dinâmica da Pressão P(z,r,t)
for p=1:pf-1
figure(6)

subplot(1,1,1)
plot(z,P(:,p)/100000)
title('Pressão no reator','FontSize',14)
xlabel('Z comprimento (m)')
ylabel('P (bar)')
xlim([0 L])
ylim([Pmin Pmax]/100000)
drawnow
end

%% Gráfico do Eta local
figure(7)
plot(rp,eta_local(2,:))
title('Fator de Efetividade Interna','FontSize',14)
xlabel('R raio do Pellet (m)')
ylabel('Eta')
xlim([0 Rp])
ylim([0 1])
drawnow
% surf(rp,z,eta)     


