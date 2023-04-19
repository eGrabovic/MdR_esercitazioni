%% esercitazione dinamica seriali: computed torque control
clc; clear all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

addpath(genpath('../Casadi'));
addpath(genpath('Data'));
addpath(genpath('utils'));

%% redundant robot parametrization

% robot KUKA 9R

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

% build functions from symbolic expressions
T0E_fun = Function('T0E', {q}, {T0E});
Tj_fun = Function('Tj', {q}, {Tj{:}});
Jac_fun = Function('Jac', {q}, {Jac});


%% calcolo matrici della dinamica forma standard
[B, C, G, B_dq] = stdDynFromDH(DH_table, Jtype_list, q, qdot, cg_list, mass_list, inertia_list, 'baseOffset', Toffset0, 'eeOffset', ToffsetE, 'gravity', [0;0;-9.81]);

% volendo potrei calcolare la dinamica di un sistema con incertezze sulle
% inerzie per vedere come si comporta la legge computed torque
% [B_real, C_real, G_real, B_dq_real] = stdDynFromDH(DH_table, Jtype_list, q, qdot, cg_list_real, mass_list_real, inertia_list_real, 'baseOffset', Toffset0, 'eeOffset', ToffsetE, 'gravity', [0;0;-9.81]);

%% caricamento dei risultati della cinematica CLIK
load('q_CLIK.mat'); % caricato nella variabile "qsol"
t = qsol(end, :); % tempo
qsol = qsol(1:end-1,:); % togliamo il tempo

import casadi.*
t_sym = MX.sym('t');
q_des = MX(nj, 1);
q_d_des = MX(nj, 1);
q_dd_des = MX(nj, 1);

% interpolazione con spline casadi

for i = 1:nj
    q_interp = casadi.interpolant('interp', 'bspline', {t}, qsol(i,:));
    q_des(i) = q_interp(t_sym);
    q_d_des(i) = jacobian(q_des(i), t_sym);
    q_dd_des(i) = jacobian(q_d_des(i), t_sym);
end

q_des_fun = Function('qd', {t_sym}, {q_des});
q_d_des_fun = Function('qd_d', {t_sym}, {q_d_des});
q_dd_des_fun = Function('qd_dd', {t_sym}, {q_dd_des});

%% simulazione dinamica
% stati sistema
s = [q; qdot];

% ingressi sistema
u = casadi.MX.sym('u', 6, 1);

% computed torque
omega_n = 10;
Kp = eye(nj).*omega_n.^2;
Kv = eye(nj).*2.*omega_n;
tau = B*(q_dd_des - Kv*(qdot - q_d_des) - Kp*(q - q_des)) + C*qdot + G;

tau_disturbed = B*(q_dd_des - Kv*(qdot - q_d_des) - Kp*(q - q_des)) + C*qdot + G + Jac'*u;

tau_fun = casadi.Function('tau', {t_sym, q, qdot}, {tau});
tau_fun = casadi.Function('tau', {t_sym, q, qdot}, {jacobian(tau, t_sym)});
s_dot = [qdot; inv(B)*( -C*qdot - G + tau_disturbed)];

% parametri integrazione
t0 = t(1);
tf = t(end);
dt = (tf-t0)./(length(t)-1);
Nd = length(t0:dt:tf);

% condizioni iniziali
q0 = full(q_des_fun(0));
q0_dot = zeros(nj, 1); %full(q_d_des_fun(0)).*0;
s0 = [q0; q0_dot];

% integrazione con Runge-Kutta
% aggiungo una perturbazione sull end-effector
u_num = zeros(6, Nd);
u_num(:, floor(Nd/2):floor(Nd/2)+20) = repmat([0;0;50;0;0;0], 1, 21);
sol = RK4(s, u, s_dot, dt, tf, s0, t0, u_num, 't_expr', t_sym, 'NstepsInterval', 2);

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

set(gca, 'zlim', [-1.5 1]);
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
xlabel('t s')
ylabel('$q_\textrm{i}$')

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