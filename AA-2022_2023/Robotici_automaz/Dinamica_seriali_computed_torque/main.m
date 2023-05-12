%% esercitazione dinamica di un robot seriale
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

addpath(genpath('../Data'));
addpath(genpath('utils'));


% definizione robot
nj = 9; % number of joints

dist = 0.2;
a = 0.05;

%           a  alpha      d        theta
DH_table = [0, -pi/2,    dist,            0;...
            dist,  0,       0,        -pi/2;...
            dist,  0,       0,            0;...
            dist,  0,       0,            0;...
            dist,  0,       0,            0;...
            dist,  0,       0,            0;...
            dist,  0,       0,            0;...
            dist,  0,       0,            0;...
            dist,  0,       0,            0];

DH_table = DH_table(1:nj, :);

Jtype_list(1:nj) = 'R';
b_list = -[dist, 0, 0, 0, 0, 0, 0, 0, 0, 0]./2; % joint position w.r.t. DH frame


% definizione proprietà inerziali
mass_list = 3.*ones(1, nj); % masse dei link

inertia_list = cell(1, nj);
for j = 1:nj
    if j == 1
        inertia_list{j} = mass_list(j).*diag([1/12.*dist.^2., 1/6.*a.^2, 1/12.*dist.^2]); % tensori di inerzia
    else
        inertia_list{j} = mass_list(j).*diag([1/6.*a.^2, 1/12.*dist.^2, 1/12.*dist.^2]); % tensori di inerzia
    end
end

cg_list = [[0;dist./2;0], -[dist./2;0;0],...
          -[dist./2;0;0], -[dist./2;0;0],...
          -[dist./2;0;0], -[dist./2;0;0],...
          -[dist./2;0;0], -[dist./2;0;0],...
          -[dist./2;0;0]];

%% calcolo cinematica
import casadi.*
q = MX.sym('q', nj, 1);  % symbolic variables initialization
qdot = MX.sym('qdot', nj, 1);

% symbolic forward kinematics with DH template transforms
[T0E, Tj, T0j] = DHFWkin(DH_table, q, Jtype_list);

% base offset, if any
Toffset0 = Ttz(1)*TrotX(pi);
ToffsetE = TrotY(pi/2);
Jac = DHJac(T0j, Jtype_list, 'baseOffset', Toffset0, 'eeOffset', ToffsetE);

% EE offset and base offset, if any
T0E = Toffset0*T0E;
T0E = T0E*ToffsetE;

% build functions from symbolic expressions
T0E_fun = Function('pippo', {q}, {T0E});
Tj_fun = Function('Tj', {q}, {Tj{:}});
Jac_fun = Function('Jac', {q}, {Jac});

%% calcolo matrici della dinamica forma standard
[B, C, G, B_dq] = stdDynFromDH(DH_table, Jtype_list, q, qdot, cg_list, mass_list, inertia_list, 'baseOffset', Toffset0, 'eeOffset', ToffsetE, 'gravity', [0;0;-9.81]);

%% simulazione dinamica
% stati sistema
s = [q; qdot];

% dinamica sistema
wEE = [200;0;0;0;0;0]; % wrench end-effector
tau = -0*linspace(2, 0.6, nj)'.*qdot; % smorzamento
tau = tau;% + Jac'*wEE;

q_dotdot = inv(B)*(-C*qdot - G + tau);
s_dot = [qdot; q_dotdot];

% parametri integrazione
t0 = 0;
tf = 20;
dt = 0.005;

% condizioni iniziali
q0 = zeros(nj, 1) + pi/4;
q0_dot = zeros(nj, 1);
s0 = [q0; q0_dot];

% integrazione con Runge-Kutta
sol = RK4(s, [], s_dot, dt, tf, s0, t0, []);

% example with built-in matlab integrators (remember to use t_num of this output!)
% s_dot_fun = casadi.Function('sdot', {s}, {s_dot});
% [t_num, sol] = ode15s(@(t, x) full(s_dot_fun(x)), [t0 tf], s0);
% sol = sol';

%% grafica
T0E_num = full(T0E_fun(q0));
Tj_num = cell(nj, 1);
[Tj_num{:}] = Tj_fun(q0);
Tj_num = cellfun(@full, Tj_num, 'UniformOutput',false);

transforms = init_plot_DH(Tj_num, b_list, Jtype_list,  'joint_len', 0.06, 'joint_r', 0.02, 'T0', Toffset0, 'TE', ToffsetE);

% scene plot options
view(-35,30)
xlabel('x')
ylabel('y')
zlabel('z')

set(gca, 'zlim', [-0.5 1.5]);
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

%% plot risultati
qsol = full(sol(1:nj, :));
qsol_dot = full(sol(nj+1:nj*2, :));
t_num = t0:dt:tf;
figure; hold on

%plot configurazione giunti
ax = subplot(2,1,1);
plot(ax, t_num, qsol)
xlabel('t (s)')
ylabel('$q_i$')
set(gca, 'Fontsize', 26)

%plot velocità giunti
ax = subplot(2,1,2);
plot(ax, t_num, qsol_dot)
xlabel('t (s)')
ylabel('$\dot{q}_i$')
set(gca, 'Fontsize', 26)

%% verifica teorema delle forze vive
% dT/dt - (tau - G)qdot =0

% derivata energia cinetica
Bdot = zeros(nj, nj, 'casadi.MX');
for j = 1:nj
    Bdot = Bdot + B_dq{j}.*qdot(j);
end
dT_dt_fun = casadi.Function('T', {q, qdot}, {qdot'*B*q_dotdot + 0.5*qdot'*Bdot*qdot});
dT_dt_num = full(dT_dt_fun(qsol, qsol_dot));

% potenza immessa nel sistema
P_fun = casadi.Function('P', {q, qdot}, {(tau - G)'*qdot});
P_num = full(P_fun(qsol, qsol_dot));

diff = dT_dt_num - P_num;

figure; hold on;
plot(t_num, diff)
fprintf('max difference %f\n', max(abs(diff)));
set(gca, 'Fontsize', 26)

%% calcolo energia meccanica

% en. cinetica
T_fun = casadi.Function('T', {q, qdot}, {0.5*qdot'*B*qdot});
T_num = full(T_fun(qsol, qsol_dot));

% en. potenziale
U_expr = 0;
for j = 1:nj
    pcj = Toffset0*T0j{j}*[cg_list(:, j);1];
    U_expr = U_expr - mass_list(j).*[0, 0, -9.81]*pcj(1:3);
end
U_fun = casadi.Function('U', {q}, {U_expr});
U_num = full(U_fun(qsol));

% en. mecc 
EM = T_num + U_num;
figure; hold on
plot(t_num, T_num, 'linewidth', 1.4)
plot(t_num, U_num, 'linewidth', 1.4)
plot(t_num, EM, 'linewidth', 1.4)

xlabel('t (s)')
ylabel('(J)')
legend('$T(t)$', '$U(t)$', '$E_\textrm{M}(t)$')
set(gca, 'fontsize', 26)
