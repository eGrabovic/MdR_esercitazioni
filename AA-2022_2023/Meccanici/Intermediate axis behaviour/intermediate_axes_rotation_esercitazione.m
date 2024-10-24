clc; clear all;

%% comandi preliminari
% comandi per cambiare la formattazione grafica delle scritte nei grafici
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% aggiungo al path tutte le cartelle necessarie
addpath(genpath(pwd))
%addpath(genpath('..\Casadi'))
addpath('C:\casadiPackages\casadi-windows-matlabR2016a-v3.5.5\')

%% definizione delle proprietà del corpo

% tensore d'inerzia
J = eye(3);
% impostiamo il tensore d'inerzia in modo che l'asse di rotazione z abbia
% il momento intermedio
J(1, 1) = 0.002*J(1,1);  % inerzia asse X (inerzia minore)
J(2, 2) = 0.009*J(2,2);  % inerzia asse Y (inerzia maggiore)
J(3, 3) = 0.006*J(3,3);  % inerzia asse Z (asse intermedio)
J(1, 3) = 1e-5;   % piccola perturbazione al tensore di inerzia
J(3, 1) = 1e-5;

% postura iniziale (orientazione)
q0 = [1;0;0;0];
% velocità angolare iniziale
omega0 = [0;0.01;2];

% omega0 = [20;0;0];

%% scriviamo la dinamica con CasADi
import casadi.*
q = SX.sym('q', 4, 1); % parametri di Eulero
omega = SX.sym('omega', 3, 1); % velocità angolare
u = SX.sym('u', 3, 1); % momenti esterni

% stati
x = vertcat(q, omega);

% equazioni che collegano la velocità angolare con la derivata dei
% parametri di Eulero
% 
stab = 1;

qdot = EulParBodyJacInv(q)*omega + stab.*(1./sum(q.^2) - 1).*q; % stabilizzazione di Baumgarte per mantenere q nel manifold dei quaternioni unitari

% scriviamo le equazioni di Eulero
omegadot =  (J)\(u - hat(omega)*(J)*omega);

%omegadot =  inv(J)(u - hat(omega)*(J)*omega);

% dinamica nello spazio degli stati
xdot = [qdot; omegadot];

%% Analisi della stabilizzazione di Baumgarte

qtest = [0.5,0.6,0.8,1]; % quaternione con norma ||q|| = 1.5
% guardiamo cosa fa il termine di stabilizzazione a tale quaternione dopo
% N = 15 passi
N = 15;
qti = nan(1,N);
qti(1) = norm(qtest);
for s = 1:N
    qtest = ( + stab.*(1./sum(qtest.^2) - 1).*qtest).*0.05 + qtest; 
    qti(s+1) = norm(qtest);
end

figure('color', 'w'); hold on; grid on; box on;

line(1:N+1, qti, 'color', 'k', 'linewidth', 1.4)
xlabel('t')
ylabel('$||q||$')
set(gca, 'Fontsize', 30)

%% Simulazione (integrazione numerica)
dt = 0.05; % step di integrazione
t_end = 60; % tempo finale

u_num = zeros(3, t_end./dt+1);
% coppia di perturbazione
% u_num(:, round(end/3): round(end/3)+3) = repmat([0; 0.04; 0], 1, 4);

% integriamo con schema Runge-Kutta
x_num = RK4(x, u, xdot, dt, t_end, [q0; omega0], 0, u_num);

% integriamo con Eulero in avanti
% x_num = forwardEuler(x, u, xdot, dt, t_end, [q0; omega0], 0, u_num);

% integriamo con Eulero indietro
% x_num = backwardsEuler(x, u, xdot, dt, t_end, [q0; omega0], 0, u_num);

% estraiamo i parametri di Eulero per l'animazione
q_num = x_num(1:4,:);
q_norm = (vecnorm(q_num, 2, 1)); % calcoliamo la norma dei quaternioni 
q_num = q_num./q_norm; % si toglie sporcizia numerica che non rende q unitario

%% Animazione

% inizializziamo la scena grafica
R = quatToRot(q0);
d = [0;0;0];
G = RpTohomogeneous(R, d);
[ax, transform] = init_plot_body(G);

for i = 1:length(q_num)
   R = quatToRot(q_num(:, i)); 
   transform.Matrix = RpTohomogeneous(R, d);
   pause(0.005)
    drawnow
end

%% plot dei risultati

% estraiamo la velocità angolare
omega_num = x_num(5:7,:);

t = 0:dt:t_end;
figure('color', 'w'); hold on; box on; grid on
xlabel('t')
ylabel('$\omega$')

line(t, omega_num(1,:), 'color', 'b', 'linewidth', 1.4);
line(t, omega_num(2,:), 'color', 'r', 'linewidth', 1.4);
line(t, omega_num(3,:), 'color', 'k', 'linewidth', 1.4);
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$'})
set(gca, 'Fontsize', 30)

% colcolo energia cinetica
Ec = 0.5.*(dot(omega_num, J*omega_num, 1));

figure('color', 'w'); hold on; box on; grid on
title('Energia cinetica del sistema')
xlabel('t')
ylabel('$E_c$')

line(t, Ec, 'color', 'b', 'linewidth', 1.4);
set(gca, 'Fontsize', 30)

% plot dell'errore del quaternione 
figure('color', 'w'); hold on; box on; grid on
xlabel('t')
ylabel('$||q||-1$')
line(t, q_norm-1, 'color', 'b', 'linewidth', 1.4);
set(gca, 'Fontsize', 30)

%% integrazione con la dinamica scritta in componenti spatial
import casadi.*
q = SX.sym('q', 4, 1);
omega = SX.sym('omega', 3, 1);
u = SX.sym('u', 3, 1);

% stati
x = vertcat(q, omega);

% scriviamo le equazioni che collegano la velocità angolare con la derivata dei parametri
% di Eulero
qdot = EulParSpatialJacInv(q)*omega + stab.*(1./sum(q.^2) - 1).*q;

% scriviamo le equazioni di Eulero
M = zeros(3,1); % momenti esterni
Rsb = quatToRot(q);
Js = Rsb*J*Rsb.';
omegadot =  (Js)\(M - hat(omega)*(Js)*omega);

% dinamica nello spazio degli stati
xdot = [qdot; omegadot];

dt = 0.05;
t_end = 60;

% integriamo con schema Runge-Kutta
x_num = RK4(x, [], xdot, dt, t_end, [q0; omega0], 0, []);

% integraiamo con Eulero in avanti
% x_num = forwardEuler(x, [], xdot, dt, t_end, [q0;omega0], 0, []);

% integriamo con Eulero indietro
% x_num = backwardsEuler(x, [], xdot, dt, t_end, [q0;omega0], 0, []);

% estraiamo i parametri di Eulero per l'animazione
q_num = x_num(1:4,:);

% inizializziamo la scena grafica
G = RpTohomogeneous(quatToRot(q0), d);
[ax, transform] = init_plot_body(G);

q_num = q_num./(vecnorm(q_num, 2, 1));
for i = 1:length(q_num)
   transform.Matrix = RpTohomogeneous(quatToRot(q_num(:, i)), d);
   drawnow
end

omega_num = x_num(5:7,:);

t = 0:dt:t_end;

% plot velocità angolari
figure('color', 'w'); hold on; box on; grid on
xlabel('t')
ylabel('$\omega$')

line(t, omega_num(1,:), 'color', 'b', 'linewidth', 1.4);
line(t, omega_num(2,:), 'color', 'r', 'linewidth', 1.4);
line(t, omega_num(3,:), 'color', 'k', 'linewidth', 1.4);
legend({'$\omega_x$', '$\omega_y$', '$\omega_z$'})
set(gca, 'Fontsize', 30)

% calcolo dell'energia cinetica
Ec = 1:length(q_num);
for i = 1:length(q_num)
   R = quatToRot(q_num(:, i)); 
   Ec(i) = 0.5*omega_num(:, i)'*(R*J*R')*omega_num(:, i);
end

% plot energia cinetica
figure('color', 'w'); hold on; box on; grid on
title('Energia cinetica del sistema')
xlabel('t')
ylabel('$E_c$')

line(t, Ec, 'color', 'b', 'linewidth', 1.4);
set(gca, 'Fontsize', 30)