%% esercitazione dinamica vincolata
clc; clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

%! ricordare di aggiungere CasADi al path di matlab
addpath(genpath('../utils'));
addpath(genpath('../graphics'));

%% esempio semplice di warm up: quadrilatero articolato nel piano

% il quadrilatero lo componiamo chiudendo insieme 2 RR planari:
% il primo RR viene definito con etichetta A mentre il secondo con
% etichetta B

% proprietà geometriche link
dA(1) = 0.2; % estensioni link
dA(2) = 0.1;
dB(1) = 0.15;
dB(2) = 0.3;
dAB = 0.17; % distanza tra base di A e base di B

a = 0.05; % spessori link

njA = 2;
njB = 2;

% tabelle DH
DH_tableA = [dA(1),    0,       0,        0;...
             dA(2),    0,       0,        0];

DH_tableB = [dB(1),    0,       0,        0;...
             dB(2),    0,       0,        0];

% tipologia giunti
Jtype_listA(1:njA) = 'R';
Jtype_listB(1:njB) = 'R';

% posizionamento giunto w.r.t. DH frame (solo per grafica)
b_listA = [0, 0];
b_listB = [0, 0];

% proprietà inerziali
mass_listA = 3.*ones(1, njA); % masse dei link
inertia_listA{1} = mass_listA(1).*diag([1/6.*a.^2, 1/12.*dA(1).^2, 1/12.*dA(1).^2]);
inertia_listA{2} = mass_listA(2).*diag([1/6.*a.^2, 1/12.*dA(2).^2, 1/12.*dA(2).^2]);

mass_listB = 3.*ones(1, njA); % masse dei link
inertia_listB{1} = mass_listB(1).*diag([1/6.*a.^2, 1/12.*dB(1).^2, 1/12.*dB(1).^2]);
inertia_listB{2} = mass_listB(2).*diag([1/6.*a.^2, 1/12.*dB(2).^2, 1/12.*dB(2).^2]);

% quote baricentriche dei link in componenti locali DH
cg_listA = [[-dA(1);0;0],[-dA(2);0;0]];
cg_listB = [[-dB(1);0;0],[-dB(2);0;0]];

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
Toffset0A = TrotY(pi/2)*TrotZ(pi);
Toffset0B = Tty(dAB)*TrotY(pi/2)*TrotZ(pi);

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

%% vincolo cinematico di chiusura 
% supponiamo incastro tra i due link di estremità
% ciò significa che il twist all EE di A deve essere sempre uguale al twist
% EE di B (velocità lineare ed angolare)

% matrici di vincolo non servono: perchè?
H = eye(6); % incastro

% matrici di Grasp non servono: perchè?

% matrice di vincolo Pfaffiana
A = [H*JacA, -H*JacB]; % in pratica così l'espressione A*qdot = 0 mi indica che JacA*qdotA = JacB*qdotB

% !! Attenzione, bisogna scegliere i g.d.l. del vincolo non nulli,
% altrimenti si avrà matrice singolare A*Binv*A^T
A = A([2,3,4], :); % la nostra cinemtica è piana in z-y. Scelgo quindi vy, vz, omega_x

% derivata temporale di A
q = [qA;qB];
qdot = [qdotA;qdotB];
A_dq = jacobian(A, q);

Adot = zeros(size(A), 'casadi.SX');
for i = 1:njA+njB
    A_dqi = reshape(A_dq(:, i), size(A));
    Adot = Adot + A_dqi.*qdot(i);
end

%% dinamica vincolata: simulazione con il metodo Lagrangiano aumentato
B = blkdiag(BA, BB);
C = blkdiag(CA, CB);
G = vertcat(GA, GB);

% applicchiamo uno smorzamento ai giunti ed una molla al primo giunto di A
tauA = -0.05.*qdotA - diag([0, 0])*qA;
tauB = -0.05.*qdotB;
tau = [tauA; tauB];

% precalcolo alcune matrici utili per agevolare i calcoli
Binv = pinv(B);
AT = A';
invABAT = pinv(A*Binv*AT);
Ap = Binv*AT*invABAT;

% inversa matrice di massa del sistema aumentato
M = [(eye(njA+njB) - Ap*A)*Binv,       Ap;...
                            Ap', -invABAT];

Q = [tau - C*qdot - G;-Adot*qdot];  

qdotdot = M(1:njA+njB, :)*Q;

% attenzione i moltiplicatori non fanno parte degli stati: non bisogna
% integrarli. E' un' espressione algebrica che può essere valutata
% post-integrazione
lambda = M(njA+njB+1:end, :)*Q;
lam_fun = casadi.Function('lam', {[q; qdot]}, {lambda});

% stati
s = [q; qdot];

% dinamica stati
sdot = [qdot; qdotdot];

% parametri integrazione
t0 = 0;
tf = 15;
dt = 0.005;

q0A = [pi/4; -pi/2-pi/3];
q0B = [-pi/4; pi/3];
q0dotA = zeros(njA, 1);
q0dotB = zeros(njB, 1);

% attenzione: non posso partire con condizioni iniziali arbitrarie, devo
% rispettare il vincolo di chiusura. posso imporre qA e calcolare le qB

pos_err = T0EA(2:3, 4) - T0EB(2:3, 4);
xA = T0EA(1:3, 1);
xB = T0EB(1:3, 1);
or_err = xA(2)*xB(3) - xA(3)*xB(2); % annullare componente del prodotto vettoriale

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
problem = struct('x', [qA; qB], 'f', 0, 'g', err);

% costruzione problema di ottimizzazione. Vengono calcolate in automatico
% le derivate per valutare hessiano e gradiente del Lagrangiano. Viene
% stabilita un' interfaccia col solutore IPOPT (algoritmo interior-point, stato dell'arte)
S = casadi.nlpsol('S', 'ipopt', problem);

% chiamata al solutore con initial guess numeriche e bound sui vincoli
sol = S('x0', [q0A; q0B], 'ubg', 0, 'lbg', 0); % upper e lower bound 0 sui vincoli significa che sono equality

% estrazione della soluzione
q0 = full(sol.x);
q0A = q0(1:2);
q0B = q0(3:4);

s0 = vertcat(q0, q0dotA, q0dotB);

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
qsolA = full(sol(1:njA, :));
qsolB = full(sol(njA+1:njA+njB, :));

for j = 1:size(qsolA, 2)
    
   updatePlot_DH(Tj_funA, qsolA(:, j), Toffset0A, ToffsetEA, transformsA)  % update manipulator joints graphics
   updatePlot_DH(Tj_funB, qsolB(:, j), Toffset0B, ToffsetEB, transformsB)  % update manipulator joints graphics
   pause(0)
end