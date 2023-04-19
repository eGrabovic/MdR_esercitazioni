%% esercitazione dinamica vincolata
clc; clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

addpath(genpath('../Casadi'));
addpath(genpath('../utils'));

%% esempio semplice di warm up: meccanismo biella manovella

% il quadrilatero lo componiamo chiudendo insieme 2 RR planari:
% il primo RR viene definito con etichetta A mentre il secondo con
% etichetta B

% proprietà geometriche link
d(1) = 0.1;
d(2) = 0.25;
d(3) = 0.8;

a = 0.05; % spessori link

nj = 3;

% tabelle DH
DH_table = [d(1),    0,       0,        0;...
            d(2),    0,       0,        0;...
               0,    0,       0,        0];


% tipologia giunti
Jtype_list(1:nj) = 'R';

% posizionamento giunto w.r.t. DH frame (solo per grafica)
b_list = [0, 0, 0];

% proprietà inerziali
mass_list = 2.*ones(1, nj); % masse dei link
mass_list(1) = 5;
inertia_list{1} = mass_list(1).*diag([1/6.*a.^2, 1/12.*d(1).^2, 1/12.*d(1).^2]);
inertia_list{2} = mass_list(2).*diag([1/6.*a.^2, 1/12.*d(2).^2, 1/12.*d(2).^2]);
inertia_list{3} = mass_list(3).*diag([1/6.*a.^2, 1/12.*d(3).^2, 1/12.*d(3).^2]);

% quote baricentriche dei link in componenti locali DH
cg_list = [[-d(1)*2;0;0],[-d(2)/2;0;0], [0;0;0]];

%% calcolo cinematica
import casadi.*
q = SX.sym('qA', nj, 1);  % symbolic variables initialization
qdot = SX.sym('qdot', nj, 1);

% symbolic forward kinematics with DH template transforms
[T0E, Tj, T0j] = DHFWkin(DH_table, q, Jtype_list);

% base offset, if any
Toffset0 = TrotY(pi/2)*TrotZ(pi);

% EE offset and base offset, if any
ToffsetE = TrotY(pi/2);

T0E = Toffset0*T0E*ToffsetE;

% jacobians calculation
Jac = DHJac(T0j, Jtype_list, 'baseOffset', Toffset0, 'eeOffset', ToffsetE);

% build functions from symbolic expressions
T0E_fun = Function('T0EA', {q}, {T0E});
Tj_fun = Function('TjA', {q}, {Tj{:}});

%% calcolo matrici della dinamica forma standard
[B, C, G] = stdDynFromDH(DH_table, Jtype_list, q, qdot, cg_list, mass_list, inertia_list, 'baseOffset', Toffset0, 'eeOffset', ToffsetE, 'gravity', [0;9.81;0]);

%% vincolo cinematico di chiusura 
% supponiamo incastro tra i due link di estremità
% ciò significa che il twist all EE di A deve essere sempre uguale al twist
% EE di B (velocità lineare ed angolare)

H = null([0;1;0;0;0;0]')';
% matrice di vincolo Pfaffiana
A = H*Jac;

A = A([2,3], :);% la nostra cinemtica è piana in z-y. la velocità lungo z e rotazione attorno a x devono essere nulli

% derivata temporale di A
A_dq = jacobian(A, q);

Adot = zeros(size(A), 'casadi.SX');
for i = 1:nj
    A_dqi = reshape(A_dq(:, i), size(A));
    Adot = Adot + A_dqi.*qdot(i);
end

%% dinamica vincolata: simulazione con il metodo Lagrangiano aumentato

% applicchiamo uno smorzamento ai giunti ed una molla al primo giunto di A
t_sym = casadi.SX.sym('t', 1, 1);
tau = -0.02.*qdot; % smorzamento sui giunti
delta = T0E(2,4)+d(1)+d(2); % indicatore di posizione del pistone
tau = tau + Jac'*[0; 1000.*(t_sym < 0.015) + 300.*(delta <= 0.015 & cos(q(1)) < -0.0001); 0; 0; 0; 0]; % simuliamo un'esplosione sul pistone quando è al punto superiore

% precalcolo alcune matrici utili per agevolare i calcoli
Binv = pinv(B);
AT = A';
invABAT = pinv(A*Binv*AT);
Ap = Binv*AT*invABAT;

% inversa matrice di massa del sistema aumentato
M = [(eye(nj) - Ap*A)*Binv,       Ap;...
                       Ap', -invABAT];

Q = [tau - C*qdot - G;-Adot*qdot];                        
qdotdot = M(1:nj, :)*Q;

% attenzione i moltiplicatori non fanno parte degli stati: non bisogna
% integrarli. E' un' espressione algebrica che può essere valutata
% post-integrazione
lambda = M(nj+1:end, :)*Q;
lam_fun = casadi.Function('lam', {[q; qdot]}, {lambda});

% stati
s = [q; qdot];

% dinamica stati
sdot = [qdot; qdotdot];

% parametri integrazione

t0 = 0;
tf = 30;
dt = 0.005;

q0 = [pi/2-0.05; pi/4; 0];
q0dot = zeros(nj, 1);

% attenzione: non posso partire con condizioni iniziali arbitrarie, devo
% rispettare il vincolo di chiusura. posso imporre qA e calcolare le qB

pos_err = T0E(3, 4);
x = T0E(1:3, 3);
xcon = [0;1;0];
or_err = x'*xcon-1;

err = vertcat(pos_err, or_err);

% in pratica, avendo 3 condizioni sull'errore di chiusura del vincolo
% dovrei imporre una variable di giunto come nota (i.e. qA_1) e fare un rootfinding
%(soluzione sistema non lineare) per trovare le qB e qA_2. Il problema è che anche la qA_1 non può essere
% scelta arbitrariamente (per certi suoi valori non esistono soluzioni
% degli altri giunti che rispettano il vincolo). Per ovviare a questo problema
% ed avere una inizializzazione robusta uso un algoritmo di ottimizzazione
% dove si identificano contemporaneamente tutti valori dei giunti di A e B 
% mentre si impone l'annullamento dell'errore come vincolo. Non c'è una
% funzione di costo da minimizzare: tali problemi si definiscono 'problemi
% di feasibility'. Trovano, tra le infinite soluzioni che rispettano il
% vincolo, quella che è più vicina ai valori di initial guess.

% struttura dati che definisce il problema:
% 'x' -> variabili di ottimizzazione;
% 'f' -> funzione di costo da minimizzare;
% 'g' -> vincoli (nonlineari e lineari; equality e unequality) da far rispettare
problem = struct('x', [q], 'f', 0, 'g', err);

% costruzione problema di ottimizzazione. Vengono calcolate in automatico
% le derivate per valutare hessiano e gradiente del Lagrangiano. Viene
% stabilita un' interfaccia col solutore IPOPT (algoritmo interior-point, stato dell'arte)
S = casadi.nlpsol('S', 'ipopt', problem);

% chiamata al solutore con initial guess numeriche e bound sui vincoli
sol = S('x0', [q0], 'ubg', 0, 'lbg', 0); % upper e lower bound 0 sui vincoli significa che sono equality

% estrazione della soluzione
q0 = full(sol.x);

s0 = vertcat(q0, q0dot);

sol = RK4(s, [], sdot, dt, tf, s0, t0, [], 't_expr', t_sym);

%% proviamo a risolvere con integratore a step adattivo (ode45 di Matlab)

sdot_fun = casadi.Function('sdot', {t_sym, s}, {sdot});

sol_ode45 = ode45(@(t, s) full(sdot_fun(t, s)), [t0 tf], s0);
% ricordarsi full() per convertire DM casadi in double
% oltre al ode45 ci sono altri integratori per problemi 'stiff'.
% provare con tali solver e confrontare le soluzioni

% confronto ode e RK4 a step costante
figure('color', 'w'); hold on
plot(sol_ode45.x, sol_ode45.y(1, :), 'linewidth', 1.4)
plot(t0:dt:tf, sol(1,:), 'linewidth', 1.4)
legend('ode45 (Matlab)', 'fixed step RK4')
xlabel('$t$ (s)')
ylabel('$\dot{q}_1(t)$ $(\frac{rad}{s})$')
set(gca, 'Fontsize', 26)

%% grafica
T0E_num = full(T0E_fun(q0));
Tj_num = cell(nj, 1);
[Tj_num{:}] = Tj_fun(q0);
Tj_num = cellfun(@full, Tj_num, 'UniformOutput',false);


[transforms, ax] = init_plot_DH(Tj_num, b_list, Jtype_list,...
    'joint_len', 0.06,...
    'joint_r', 0.02,...
    'T0', Toffset0,...
    'TE', ToffsetE);

plotFrame(eye(4), 'scale', 0.5, 'label', '0');
% scene plot options
view(90,0)
xlabel('x')
ylabel('y')
zlabel('z')

set(gca, 'zlim', [-1.0 1.0]);
set(gca, 'xlim', [-1.0 1.0]);
set(gca, 'ylim', [-1.0 1.0]);
set(gca, 'box', 'on')

%% animazione
qsol = full(sol(1:nj, :));
qsoldot = full(sol(nj+1:nj*2, :));
t = full(sol(end, :));

% scommentare qui per animazione con soluzione dell ode45
% qsol = sol_ode45.y(1:nj, :);

for j = 1:size(qsol, 2)
    
   updatePlot_DH(Tj_fun, qsol(:, j), Toffset0, ToffsetE, transforms)  % update manipulator joints graphics
   pause(0)
   
end