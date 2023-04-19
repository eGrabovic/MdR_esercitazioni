%% esercitazione dinamica vincolata
clc; clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

addpath(genpath('../Casadi'));
addpath(genpath('../utils'));

%% manipolazione di corpo rigido con seriali cooperanti

% abbiamo il seriale con etichetta A ed il seriale con etichetta B che
% cooperano per manipolare il corpo rigido R

njA = 3;
njB = 3;

% proprietà geometriche link
dA = 0.3.*ones(njA, 1); % estensioni link
dB = 0.3.*ones(njB, 1); % estensioni link

dAB = 0.17; % distanza tra base di A e base di B

a = 0.05; % spessori link



% tabelle DH
DH_tableA = [dA(1),    0,       0,        0;...
             dA(2),    0,       0,        0;...
             dA(3),    0,       0,        0];

DH_tableB = [dB(1),    0,       0,        0;...
             dB(2),    0,       0,        0;...
             dB(3),    0,       0,        0];

% tipologia giunti
Jtype_listA(1:njA) = 'R';
Jtype_listB(1:njB) = 'R';

% posizionamento giunto w.r.t. DH frame (solo per grafica)
b_listA = zeros(njA, 1);
b_listB = zeros(njB, 1);

% proprietà inerziali
mass_listA = 3.*ones(1, njA); % masse dei link
inertia_listA = cell(1, njA);
for i = 1:njA
    inertia_listA{i} = mass_listA(i).*diag([1/6.*a.^2, 1/12.*dA(i).^2, 1/12.*dA(i).^2]);
end

mass_listB = 3.*ones(1, njB); % masse dei link
inertia_listB = cell(1, njB);
for i = 1:njB
    inertia_listB{i} = mass_listB(i).*diag([1/6.*a.^2, 1/12.*dB(i).^2, 1/12.*dB(i).^2]);
end

% quote baricentriche dei link in componenti locali DH
cg_listA = [-dA'; zeros(2, njA)];
cg_listB = [-dB'; zeros(2, njB)];

%% calcolo cinematica
import casadi.*
qA = SX.sym('qA', njA, 1);  % symbolic variables initialization
qdotA = SX.sym('qdotA', njA, 1);

qB = SX.sym('qB', njB, 1);  % symbolic variables initialization
qdotB = SX.sym('qdotB', njB, 1);

% symbolic forward kinematics with DH template transforms
[T0EA, TjA, T0jA] = DHFWkin(DH_tableA, qA, Jtype_listA);
[T0EB, TjB, T0jB] = DHFWkin(DH_tableB, qB, Jtype_listB);

% base offset, if any
Toffset0A = Tty(-dAB/2)*TrotY(pi/2)*TrotZ(pi);
Toffset0B = Tty(dAB/2)*TrotY(pi/2)*TrotZ(pi);

% EE offset and base offset, if any
ToffsetEA = eye(4);
ToffsetEB = eye(4);

T0EA = Toffset0A*T0EA*ToffsetEA;
T0EB = Toffset0B*T0EB*ToffsetEB;

% jacobians calculation
JacA = DHJac(T0jA, Jtype_listA, 'baseOffset', Toffset0A, 'eeOffset', ToffsetEA);
JacB = DHJac(T0jB, Jtype_listB, 'baseOffset', Toffset0B, 'eeOffset', ToffsetEB);

% build functions from symbolic expressions
T0E_funA = Function('T0EA', {qA}, {T0EA});
Tj_funA = Function('TjA', {qA}, {TjA{:}});
T0E_funB = Function('T0EB', {qB}, {T0EB});
Tj_funB = Function('TjB', {qB}, {TjB{:}});

%% calcolo matrici della dinamica forma standard
[BA, CA, GA] = stdDynFromDH(DH_tableA, Jtype_listA, qA, qdotA, cg_listA, mass_listA, inertia_listA, 'baseOffset', Toffset0A, 'eeOffset', ToffsetEA);

[BB, CB, GB] = stdDynFromDH(DH_tableB, Jtype_listB, qB, qdotB, cg_listB, mass_listB, inertia_listB, 'baseOffset', Toffset0B, 'eeOffset', ToffsetEB);

%% caratterizzazione cinematica e dinamica del corpo rigido

% paramaetrizzazione postura
qR = casadi.SX.sym('qR', 6, 1); % sei parametri 

% orientazione yaw-pitch-roll
R0R = rotZ(qR(4))*rotY(qR(5))*rotX(qR(6));
T0R = [R0R, qR(1:3);...
        0, 0, 0, 1];

% jacobiano parametrizzazione
% per rimanere coerenti con le parametrizzazioni dei seriali 
% si sceglie come polo l'origine di {R} e componenti spatial ({S0})
% in particolare ci servirà l'inversa dello jacobiano per le equazioni di
% congruenza tra stati
JacR_inv = [eye(3, 'casadi.SX'),     zeros(3, 3, 'casadi.SX');...
       zeros(3, 3, 'casadi.SX'), spatialJac_ZYX_inv(qR(4:6))];

% creiamo funzioni valutabili numericamente
T0R_fun = casadi.Function('T0R', {qR}, {T0R});
JacR_inv_fun = casadi.Function('JacR_inv_fun', {qR}, {JacR_inv});

% parametri geometrici
L = 0.3; % lunghezza
t = 0.05; % spessore

% parametri inerziali
massR = 3;
inertiaR = massR.*diag([1/12*L^2, 1/6*t^2, 1/12*L^2]); % tens. baricentrico princp.

% è più conveniente scegliere come stati per un c.r. direttamente le
% velocità invece delle derivate della parametrizzazione

twistR = casadi.SX.sym('twistR', 6, 1);

Fext = massR.*[0;0;-9.81]; % più altre eventuali forze esterne sul c.r.
Mext = [0;0;0]; % eventuali momenti esterni agenti sul c.r.
omega = twistR(4:6);
inertia0 = R0R*inertiaR*(R0R.');

qdotR = JacR_inv*twistR; % ricostruzione stati twist-derivate param.

BR = blkdiag(massR*eye(3), inertia0); % tensore d'inerzia generalizzato
RHSR = [Fext;... % Newton
    (Mext - hat(omega)*inertia0*omega)]; % Eulero

%% vincolo cinematico di chiusura 
% supponiamo incastro tra i due link di estremità
% ciò significa che il twist all EE di A deve essere sempre uguale al twist
% EE di B (velocità lineare ed angolare)

% lavoriamo nel piano y-z e rotazioni attorno a x

% matrici di vincolo
HA = [0,1,0,0,0,0;...
      0,0,1,0,0,0];
  
HB = HA;
H = blkdiag(HA, HB);
J = blkdiag(JacA, JacB);

% matrici di grasp
PA = T0R(1:3, 1:3)*[0;-L/2;0];
PB = T0R(1:3, 1:3)*[0;+L/2;0];
graspA = graspMatrix(PA(1:3));
graspB = graspMatrix(PB(1:3));
Grasp = [graspA, graspB];

% matrice di vincolo Pfaffiana
A = [H*J, -H*Grasp'];
% Afun = casadi.Function('A', {[qA; qB; qR]}, {A});

% derivata temporale di A
q = [qA; qB; qR];
qdot  = [qdotA; qdotB; twistR];
qdotf = [qdotA; qdotB; qdotR];
A_dq = jacobian(A, q);

Adot = zeros(size(A), 'casadi.SX');
for i = 1:length(q)
    A_dqi = reshape(A_dq(:, i), size(A));
    Adot = Adot + A_dqi.*qdotf(i);
end

%% dinamica vincolata: simulazione con il metodo Lagrangiano aumentato
B = blkdiag(BA, BB, BR);
C = blkdiag(CA, CB);
G = vertcat(GA, GB);

% applicchiamo uno smorzamento ai giunti ed una molla al primo giunto di A
tauA = -0.25.*qdotA;
tauB = -0.25.*qdotB;
tau = [tauA; tauB];

% precalcolo alcune matrici utili per agevolare i calcoli
Binv = inv(B);
AT = A';
invABAT = inv(A*Binv*AT);
Ap = Binv*AT*invABAT;

% inversa matrice di massa del sistema aumentato
M = [(eye(njA+njB+6) - Ap*A)*Binv,       Ap;...
                              Ap', -invABAT];

Q = [tau - C*qdot(1:njA+njB) - G; RHSR; -Adot*qdot];
qddlam = M*Q;
qdotdot = qddlam(1:njA+njB+6);

% attenzione i moltiplicatori non fanno parte degli stati: non bisogna
% integrarli. E' un' espressione algebrica che può essere valutata
% post-integrazione
lambda = qddlam(njA+njB+7:end);
lam_fun = casadi.Function('lam', {[q; qdot]}, {lambda});

% stati
s = [q; qdot];

% dinamica stati
sdot = [qdotA; qdotB; qdotR; qdotdot];

% parametri integrazione

t0 = 0;
tf = 15;
dt = 0.005;

q0A = zeros(njA, 1)-0.5;
q0B = zeros(njB, 1)+0.1;
q0R = [0; -0.5; 0.05; 0; 0; 0.3];
q0dotA = zeros(njA, 1);
q0dotB = zeros(njB, 1);
q0dotR = zeros(6, 1);

% attenzione: non posso partire con condizioni iniziali arbitrarie, devo
% rispettare il vincolo di chiusura.

PA0 = T0R*[0;-L/2;0;1];
PB0 = T0R*[0;+L/2;0;1];
pos_errA = T0EA(2:3, 4) - PA0(2:3);
pos_errB = T0EB(2:3, 4) - PB0(2:3);

err = vertcat(pos_errA, pos_errB);

% struttura dati che definisce il problema:
% 'x' -> variabili di ottimizzazione;
% 'f' -> funzione di costo da minimizzare;
% 'g' -> vincoli (nonlineari e lineari; equality e unequality) da far rispettare
problem = struct('x', [qA; qB; qR(1:3)], 'p', qR(4:6), 'f', 0, 'g', err);

% costruzione problema di ottimizzazione. Vengono calcolate in automatico
% le derivate per valutare hessiano e gradiente del Lagrangiano. Viene
% stabilita un' interfaccia col solutore IPOPT (algoritmo interior-point, stato dell'arte)
S = casadi.nlpsol('S', 'ipopt', problem);

% chiamata al solutore con initial guess numeriche e bound sui vincoli
sol = S('x0', [q0A; q0B; q0R(1:3)], 'p', q0R(4:6), 'ubg', 0, 'lbg', 0); % upper e lower bound 0 sui vincoli significa che sono equality

% estrazione della soluzione
q0 = full(sol.x);
q0A = q0(1:njA);
q0B = q0(njA+1:njA+njB);
q0R(1:3) = q0(njA+njB+1:end);

s0 = vertcat([q0A; q0B; q0R], q0dotA, q0dotB, q0dotR);

sol = RK4(s, [], sdot, dt, tf, s0, t0, []);

%% grafica
T0E_numA = full(T0E_funA(q0A));
Tj_numA = cell(njA, 1);
[Tj_numA{:}] = Tj_funA(q0A);
Tj_numA = cellfun(@full, Tj_numA, 'UniformOutput',false);

T0E_numB = full(T0E_funB(q0B));
Tj_numB = cell(njB, 1);
[Tj_numB{:}] = Tj_funB(q0B);
Tj_numB = cellfun(@full, Tj_numB, 'UniformOutput',false);

T0R_num = full(T0R_fun(q0R));

[transformsA, ax] = init_plot_DH(Tj_numA, b_listA, Jtype_listA,...
    'joint_len', 0.06,...
    'joint_r', 0.02,...
    'T0', Toffset0A,...
    'TE', ToffsetEA,...
    'label', 'A');

transformsB = init_plot_DH(Tj_numB, b_listB, Jtype_listB,...
    'joint_len', 0.06,...
    'joint_r', 0.02,...
    'T0', Toffset0B,...
    'TE', ToffsetEB,...
    'label', 'B',...
    'parent', ax);

% rigid body plot
transformR = hgtransform('parent', ax);
transformR.Matrix = T0R_num;
[pointsR, facesR] = createParallelepiped(t, [L/2, L/2], [0;1;0], 0);
patch('faces', facesR, 'vertices', pointsR.',  'facecolor', 'green', 'Parent', transformR);
plotFrame(eye(4), 'scale', 0.1, 'label', 'R', 'parent', transformR);

% global frame
plotFrame(eye(4), 'scale', 0.5, 'label', '0', 'parent', ax);

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

qsolA = full(sol(1:njA, :));
qsolB = full(sol(njA+1:njA+njB, :));
qsolR = full(sol(njA+njB+1:njA+njB+6, :));
qdotsol = full(sol(njA+njB+7:end, :));
qsolRdot = qdotsol(njA+njB+1:njA+njB+6, :);

for j = 1:size(qsolA, 2)
    
   updatePlot_DH(Tj_funA, qsolA(:, j), Toffset0A, ToffsetEA, transformsA)  % update manipulator joints graphics
   updatePlot_DH(Tj_funB, qsolB(:, j), Toffset0B, ToffsetEB, transformsB)  % update manipulator joints graphics
   transformR.Matrix = full(T0R_fun(qsolR(:, j)));
   pause(0)
end