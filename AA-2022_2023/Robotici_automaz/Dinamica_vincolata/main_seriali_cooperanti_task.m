%% esercitazione dinamica vincolata
clc; clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

%! ricordare di aggiungere CasADi al path di matlab
addpath(genpath('../utils'));
addpath(genpath('../graphics'));

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

%% parametrizzazione cinematica
import casadi.*
qA = MX.sym('qA', njA, 1);  % symbolic variables initialization
qdotA = MX.sym('qdotA', njA, 1);

qB = MX.sym('qB', njB, 1);  % symbolic variables initialization
qdotB = MX.sym('qdotB', njB, 1);

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
qR = casadi.MX.sym('qR', 6, 1); % sei parametri 

% orientazione yaw-pitch-roll
R0R = rotZ(qR(4))*rotY(qR(5))*rotX(qR(6));
T0R = [R0R, qR(1:3);...
        0, 0, 0, 1];

% jacobiano parametrizzazione
% per rimanere coerenti con le parametrizzazioni dei seriali 
% si parametrizza il twist con polo l'origine di {R} e componenti spatial ({S0})
% in particolare ci servirà l'inversa dello jacobiano per le equazioni di
% congruenza tra stati
cl = class(qR);
JacR_inv = [eye(3, cl),     zeros(3, 3, cl);...
       zeros(3, 3, cl), spatialJac_ZYX_inv(qR(4:6))];

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

twistR = casadi.MX.sym('twistR', 6, 1);

Fext = massR.*[0;0;-9.81]; % più altre eventuali forze esterne sul c.r.
Mext = [0;0;0]; % eventuali momenti esterni agenti sul c.r.
omega = twistR(4:6);
inertia0 = R0R*inertiaR*(R0R.');

qdotR = JacR_inv*twistR; % ricostruzione stati twist-derivate param.

BR = blkdiag(massR*eye(3), inertia0); % tensore d'inerzia generalizzato

RHSR = [Fext;...         % Newton
    (Mext - hat(omega)*inertia0*omega)]; % Eulero

%% vincolo cinematico di chiusura 

% i due seriali sono collegati al coupler con 

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
q = [qA; qB; qR];
% Afun = casadi.Function('A', {[qA; qB; qR]}, {A});

% derivata temporale di A
qdot = [qdotA; qdotB; twistR]; % stati relative alle velocità
qdotf = [qdotA; qdotB; qdotR]; % derivate delle variabili di giunto (nota che qdotR è funzione del twistR)
A_dq = jacobian(A, q);

Adot = zeros(size(A), cl);
for i = 1:length(q)
    A_dqi = reshape(A_dq(:, i), size(A));
    Adot = Adot + A_dqi.*qdotf(i); % qui si usa qdotf
end

%% calcolo della cinematica differenziale inversa con CLIK per valutare le leggi orarie necessarie sui giunti

clik_seriali_cooperanti;

%% calcolo delle coppie ai giunti necessari per equilibrare il wrench inerziale del c.r.
% creo funzioni per rivalutare con espressioni simboliche dei giunti dipendenti dal tempo
RHSR_fun = casadi.Function('RHS', {qR, twistR}, {RHSR});
BR_fun = casadi.Function('RHS', {qR}, {BR});
Jac_fun = casadi.Function('Jac', {q}, {J});

II = BR_fun(qR_t)*twistR_dot_task - RHSR_fun(qR_t, twistR_task); % wrench inerziale 

Grasph = Grasp_fun(qR_t_fun(t_sym))*H';
wEAB = pinv(Grasph(2:4,:))*II(2:4); % soluzione norma minima
% !! attenzione in wEAB si sono rimosse dalla matrice Grasph le componenti
% fuori piano a causa di singolarità numeriche. Anche dal
% wrench inerziale II si sono tolte, per compatibilità dimensionale, le
% componenti fuori paino.

Jac_t = Jac_fun(q_des_fun(t_sym));
tau = (H*Jac_fun(q_des_fun(t_sym)))'*wEAB;
tau_fun = casadi.Function('tau_task', {t_sym}, {tau});

% problema con le spline che non valutano nel punto finale
tau_num = full(tau_fun(t_num));

figure('color', 'w'); hold on;
for i = 1:size(tau_num, 1)
   plot(t_num, tau_num(i, :), 'linewidth', 1.4);
end
ylabel('$\tau$ (Nm)');
xlabel('$t$ (s)');
legend({'$\tau_{1A}$','$\tau_{2A}$','$\tau_{3A}$','$\tau_{1B}$','$\tau_{2B}$','$\tau_{3B}$'})
set(gca, 'Fontsize', 26)

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

qsolA = full(qsol(1:njA, :));
qsolB = full(qsol(njA+1:njA+njB, :));
t_num = linspace(t0, tf, length(qsolA));

for j = 1:size(qsolA, 2)
    
   updatePlot_DH(Tj_funA, qsolA(:, j), Toffset0A, ToffsetEA, transformsA)  % update manipulator joints graphics
   updatePlot_DH(Tj_funB, qsolB(:, j), Toffset0B, ToffsetEB, transformsB)  % update manipulator joints graphics
   transformR.Matrix = full(T0R_task_fun(t_num(j)));
   pause(0)
end