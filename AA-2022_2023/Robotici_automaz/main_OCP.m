%% esercitazione dinamica seriali: controllo ottimo
clc; clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

addpath(genpath('Data'));
addpath(genpath('utils'));

%% redundant robot parametrization

% robot 9R

nj = 9; % number of joints
dist = 0.2;
a = 0.05;
%           a  alpha      d        theta

DH_table = [0, -pi/2,     dist,        0;...
            0, +pi/2,        0,        0;...
            0, -pi/2,     dist,        0;...
            0, +pi/2,        0,        0;...
            0, -pi/2,     dist,        0;...
            0, +pi/2,        0,        0;...
            0, -pi/2,     dist,        0;...
            0, +pi/2,        0,        0;...
            0,     0,     dist,        0];

Jtype_list(1:nj) = 'R';
b_list = -[dist, 0, dist, 0, dist, 0, dist, 0, dist, 0]./2; % joint position w.r.t. DH frame (for graphic purposes only)

%% definizione proprietà inerziali
mass_list = 3.*ones(1, nj); % masse dei link

inertia_list = cell(1, nj);
for j = 1:nj
    if mod(j, 2)
        inertia_list{j} = mass_list(j).*diag([1/12.*dist.^2., 1/6.*a.^2, 1/12.*dist.^2]); % tensori di inerzia
    else
        inertia_list{j} = mass_list(j).*diag([1/12.*dist.^2, 1/12.*dist.^2, 1/6.*a.^2]); % tensori di inerzia
    end
end

cg_list = [
           [0;dist./2;0], [0;0;dist./2],...
           [0;dist./2;0], [0;0;dist./2],...
           [0;dist./2;0], [0;0;dist./2],...
           [0;dist./2;0], [0;0;dist./2],...
           [0;dist./2;0]
           ];

cg_list = cg_list(:, 1:nj); % centro di massa w.r.t. to Oj (in {Sj})

%% calcolo cinematica
import casadi.*
q = MX.sym('q', nj, 1);  % symbolic variables initialization
qdot = MX.sym('qdot', nj, 1);

% symbolic forward kinematics with DH template transforms
[T0E, Tj, T0j] = DHFWkin(DH_table, q, Jtype_list);

% base offset, if any
Toffset0 = Ttz(0.8)*TrotX(pi);
ToffsetE = eye(4);
Jac = DHJac(T0j, Jtype_list, 'baseOffset', Toffset0, 'eeOffset', ToffsetE);

% EE offset and base offset, if any
T0E = Toffset0*T0E;
T0E = T0E*ToffsetE;
cg_ee = rigidInverse(ToffsetE)*[cg_list(:, nj); 1];
cg_list(:, nj) = cg_ee(1:3);

% build functions from symbolic expressions
T0E_fun = Function('T0E', {q}, {T0E});
Tj_fun = Function('Tj', {q}, {Tj{:}});
Jac_fun = Function('Jac', {q}, {Jac});

%% calcolo matrici della dinamica forma standard
[B, C, G, B_dq] = stdDynFromDH(DH_table, Jtype_list, q, qdot, cg_list, mass_list, inertia_list, 'baseOffset', Toffset0, 'eeOffset', ToffsetE, 'gravity', [0;0;-9.81]);

%% caricamento dei risultati della cinematica CLIK
load('q_CLIK.mat'); % caricato nella variabile "qsol"
t = qsol(end, :); % tempo
qsol = qsol(1:end-1,:); % togliamo il tempo

% numero di passi
N = length(qsol);

% calcolo delle matrici di trasformazione omogenee corrispondednti all' E-E
% esse definiranno i "traguardi" che l'E-E dovrà rispettare (path constraints)
% N.B. tali traguardi non hanno vincoli di tempo (non sono costretto a far svolgere il task in 30 s)
T0E_path = reshape(full(T0E_fun(qsol)), 4, 4, []); % matrice 4x4xN 
bata_path = nan(4, N);
p_path = nan(3, N);
for ii = 1:N
    
    [nT, thetaT] = rotToAxisAngle(T0E_path(1:3, 1:3));
    beta_num(:, ii) = axisAngleToQuat(thetaT, nT);
    d_num(:, ii) = T0E_path(1:3, 4);
    
end

%% scrittura equazioni della dinamica e passo di integrazione
% stati sistema
s = [q; qdot];
nx = nj*2;

% ingressi sistema
tau = casadi.MX.sym('tau', nj, 1);
nu = nj;

% variabili algebriche
h = casadi.MX.sym('h'); % passo di integrazione incognito (da determinare mediane l'algoritmo di ottimizzazione)

beta_path = casadi.MX.sym('R', 4, 1);
d_path = casadi.MX.sym('d', 3, 1);

z = vertcat(h);
nz = length(z);

% dinamica
q_dotdot = cse(inv(B)*(-C*qdot - G + tau));

% forma di stato
s_dot = cse([qdot; q_dotdot]); % cse ("common sub-expression elimination")
s_dot_fun = casadi.Function('sdot', {s, tau}, {s_dot});


% integrazione della dinamica (1 passo di integrazione) con Runge-Kutta
s_next = (RK4_step(s, tau, s_dot_fun, h));

% vincoli nonlineari (path constraints)
T0E = T0E_fun(q);
R0E = T0E(1:3, 1:3);
d0E = T0E(1:3, 4);
[nT, thetaT] = rotToAxisAngle(R0E(1:3, 1:3));
beta_0E = axisAngleToQuat(thetaT, nT);
g_path = vertcat(beta_0E(:), d0E) - vertcat(beta_path(:), d_path); % terna E-E sul path


% costruisco funzioni
s_next_fun = casadi.Function('dyn_int', {s, tau, h}, {s_next});
g_path_fun = casadi.Function('g', {s, beta_path, d_path}, {cse(g_path)});

ng = nj*2+7;

%% costruzione del problema (discretizzazione dinamica sull'orizzonte temporale)

q_guess = qsol(:, 1); % inizializzo coi risultati della CLICK
q_dot_guess = q_guess.*0;

q_lb = -ones(nj, 1).*pi;
q_ub = +ones(nj, 1).*pi;

q_dot_lb = -pi/6*ones(nj, 1); % rad/s
q_dot_ub = +pi/6*ones(nj, 1); % rad/s

S_lb = [q_lb; q_dot_lb];
S_ub = [q_ub; q_dot_ub];

S_guess = [q_guess; q_dot_guess];
S_0 = S_guess;
S_1 = S_0; % passo precedente


tau_ub = linspace(200, 100, nj)';
tau_lb = -tau_ub;
tau_guess = zeros(nj, 1);

z_guess = 1e-2;
z_lb = 1e-6;
z_ub = 1;


X = cell(1, N-1);
X_guess = zeros(nx+nu+nz, N-1);
X_lb = zeros(nx+nu+nz, N-1);
X_ub = zeros(nx+nu+nz, N-1);

g = cell(1, N-1);
g_lb = zeros(ng, N-1);
g_ub = zeros(ng, N-1);

for k = 1:N-1
    
    S = casadi.MX.sym(['S_' num2str(k)], nx, 1);
    U = casadi.MX.sym(['U_' num2str(k)], nu, 1);
    Z = casadi.MX.sym(['Z_' num2str(k)], nz, 1);
    
    X{k} = [S; U; Z];
    
    % vincolo dinamica
    g_dyn = S - s_next_fun(S_1, U, Z(1));
    
    % vincolo path
    g_path = g_path_fun(S, beta_num(:, k+1), d_num(:, k+1));
    
    % vincolo k-esimo passo di integrazione
    g{k}  = [g_dyn; g_path];
    
    % inital guess e bound
    X_guess(:, k) = [S_guess; tau_guess; z_guess];
    X_lb(:, k) = [S_lb; tau_lb; z_lb];
    X_ub(:, k) = [S_ub; tau_ub; z_ub];
    S_1 = S; % passo precedente
end

cost = Z'*Z;

%% costruzione problema ed interfacciamento con IPOPT
prob = struct('x', vertcat(X{:}), 'f', cost, 'g', vertcat(g{:}));
          
opts = struct;
opts.ipopt.max_iter = 5000;
opts.ipopt.print_level = 5;

solver = casadi.nlpsol('solver', 'ipopt', prob, opts);

%% lancio ottimizzazione
sol = solver('x0', X_guess(:), 'lbx', X_lb(:), 'ubx', X_ub(:), 'lbg', 0, 'ubg', 0);

%% grafica
T0E_num = full(T0E_fun(q0));
Tj_num = cell(nj, 1);
[Tj_num{:}] = Tj_fun(q0);
Tj_num = cellfun(@full, Tj_num, 'UniformOutput',false);

transforms = init_plot_DH(Tj_num, b_list, Jtype_list,  'joint_len', 0.06, 'joint_r', 0.02, 'T0', Toffset0, 'TE', ToffsetE);

% plot fixed frame
plotFrame(eye(4), 'label', '0', 'scale', 0.7) % plot fixed frame

% ground graphics
plane_pt = [[1;1], [1;-1], [-1;-1], [-1;1]].*[0.5;0.7];
patch(gca, plane_pt(1,:), plane_pt(2,:), zeros(1, 4), 'facecolor', [.7 .7 .7]);

% ceiling graphics
plane_pt = [[1;1], [1;-1], [-1;-1], [-1;1]].*[0.1;0.12];
patch(gca, 'XData', plane_pt(1,:), 'YData', plane_pt(2,:), 'ZData', zeros(1,4)+0.8, 'facecolor', [.7 .7 .7]);

% writing board graphics
boardOrigin = [0.35;0;0.5];
board_pt = Ttx(boardOrigin(1))*Tty(boardOrigin(2))*Ttz(boardOrigin(3))*TrotY(pi/2)*([[0.3;0.3;0;1], [0.3;-0.3;0;1], [-0.3;-0.3;0;1], [-0.3;0.3;0;1]].*[1;2;1;1]);
fill3(gca, board_pt(1,:), board_pt(2,:), board_pt(3, :), board_pt(3, :).*0, 'facecolor', '#454545');

% scene plot options
view(-45,45)
xlabel('x')
ylabel('y')
zlabel('z')

set(gca, 'zlim', [-0.01 1]);
set(gca, 'xlim', [-1.5 1.5]);
set(gca, 'ylim', [-1.5 1.5]);
set(gca, 'box', 'on')

%% animazione
qsol = full(sol(1:nj, :));
tracked_line = line(nan(1, size(qsol, 2)),nan(1, size(qsol, 2)),nan(1, size(qsol, 2)), 'color', 'g', 'linewidth', 1.4);
for j = 1:size(qsol, 2)
    
   updatePlot_DH(Tj_fun, qsol(1:nj, j), Toffset0, ToffsetE, transforms)  % update manipulator joints graphics
   T0E_num = full(T0E_fun(qsol(1:nj, j)));           % compute E-E matrix

%    append numerical values for the actual E-E trajectory
   tracked_line.XData(j) = T0E_num(1,4);
   tracked_line.YData(j) = T0E_num(2,4);
   tracked_line.ZData(j) = T0E_num(3,4);
   drawnow
end

%% plot dei risultati
qsol = full(sol(1:nj, :));
% errori di configurazione
t_num = full(sol(end, :));
q_err = full(sol(1:nj, :) - q_des_fun(t_num));
figure; hold on;
labels = cell(nj, 1);

for i = 1:nj
    plot(t, q_err(i,:), 'linewidth', 1.4);
    labels{i} = sprintf('$q_%i$', i);
end
legend(labels{:})
set(gca, 'fontsize', 26)
xlabel('t (s)')
ylabel('$q_\textrm{errors}$')

%% coppie di attuazione
qsol_dot = full(sol(nj+1:nj*2, :));
tau_num = full(tau_fun(t_num, qsol, qsol_dot));

figure; hold on;
labels = cell(nj, 1);
for i = 1:nj
    plot(t, tau_num(i,:), 'linewidth', 1.4);
    labels{i} = sprintf('$q_%i$', i);
end
legend(labels{:})
set(gca, 'fontsize', 26)
xlabel('t s')
ylabel('$\tau_\textrm{i}$')