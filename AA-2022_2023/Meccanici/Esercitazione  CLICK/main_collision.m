%% Esercitazione CLICK 1P8R con task secondario di evitare collisioni
% 
% Analisi cinematica differenziale inversa con inseguimento di traiettorie
% con algoritmo CLIK (Closed-Loop Inverse Kinematics)

clc; clear all;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

%% import coordinates for the writing task
load('coordinateX_scritta.mat'); % variabile X
load('coordinateY_scritta.mat'); % variabile Y

m = 128;  % number of DFT coefficients

% discrete Fourier transform
pk = DFT_vec([X;Y], m);

% inverse Fourier transform (the actual analytic expression of the trajectory)
p = @(t) IFST(t, pk, m);

norm_pk = vecnorm(pk, 2, 1);
[norm_pk, idorder] = sort(norm_pk, 'descend');

Nd = 2000; % number of sampling points

%% Task trajectory calculation 
import casadi.*
boardOrigin = [0.35;0;0.5];        % center of the writing board
t = SX.sym('t');                    % symbolic variable for time
ti = 0;                             % initial time
tf = 30;                            % final time
alpha = (t-ti)./(tf-ti).*2*pi;      % scaling of the Fourier transofrm period
pt_expr = IFST(alpha, pk, m, idorder).*0.0015;
pt_expr = Ttransl(boardOrigin)*TrotY(-pi/2)*TrotZ(3*pi/2)*vertcat(pt_expr, 0, 1);
pt_fun = Function('pt', {t}, {pt_expr});
writing_pts = nan(4, Nd);
t_num = linspace(ti, tf, Nd);
for i = 1:Nd
    writing_pts(:, i) = full(pt_fun(t_num(i)));
end

%% calcolo del twist
import casadi.*
% velocità lineare 
v_expr = jacobian(pt_expr, t);

orientation_type = 1;

if orientation_type == 1
    % % % ORIENTAZIONE OPZIONE 1 % % %
    % l'orientazione la vincoliamo in modo che la terna solidale al punto di
    % scrittura abbia l'asse x tangente alla curva
    
    t_expr = v_expr(1:3)./norm(v_expr(1:3)); % versore tangente
    w_expr = [-1;0;0];                       % versore binormale
    n_expr = -cross(w_expr, t_expr);         % versore normale alla curva piana
    
    
    R_expr = [t_expr, w_expr, n_expr]*rotX(pi/2);
    
elseif orientation_type == 2
    % % % ORIENTAZIONE OPZIONE 2 % % %
    % l'orientazione la vincoliamo in modo da avere pendenza della punta
    % proporzionale all'altezza del punto di scrittura rispetto al baricentro
    % (centro lavagna)
    
    beta_max = pi/6;
    beta_min = -pi/6;
    
    max_z = max(Y).*0.0015 + boardOrigin(3);
    min_z = min(Y).*0.0015 + boardOrigin(3);
    
    s = (pt_expr(3) - min_z)./(max_z - min_z);
    
    beta = (beta_max - beta_min).*s + beta_min;
    
    R_expr = rotY(-beta + pi/2);
    
end

T_expr = [R_expr, pt_expr(1:3); 0 0 0 1];

omega_expr = vecForm(  reshape(jacobian(R_expr, t), 3, 3) * (R_expr.')  );

twist_expr = vertcat(v_expr(1:3), omega_expr);

twist_fun = Function('twist', {t}, {twist_expr});

T_fun = Function('twist', {t}, {T_expr});

twist_num = full(twist_fun(linspace(ti, tf, Nd)));

% figure('color', 'w')
% ax = subplot(2,3,1);
% plot(ax, twist_num(1,:))
% xlabel('$t$ [s]');
% ylabel('$v_x$ [m/s]');
% ax = subplot(2,3,2);
% plot(ax, twist_num(2,:))
% xlabel('$t$ [s]');
% ylabel('$v_y$ [m/s]');
% ax = subplot(2,3,3);
% plot(ax, twist_num(3,:))
% xlabel('$t$ [s]');
% ylabel('$v_z$ [m/s]');
% ax = subplot(2,3,4);
% plot(ax, twist_num(4,:))
% xlabel('$t$ [s]');
% ylabel('$\omega_x$ [m/s]');
% ax = subplot(2,3,5);
% plot(ax, twist_num(5,:))
% xlabel('$t$ [s]');
% ylabel('$\omega_y$ [m/s]');
% ax = subplot(2,3,6);
% plot(ax, twist_num(6,:))
% xlabel('$t$ [s]');
% ylabel('$\omega_z$ [m/s]');
% 
% set(findall(gcf,'-property','FontSize'),'FontSize', 30);

%% redundant robot parametrization

% robot KUKA 9R

nj = 9; % number of joints

%           a  alpha      d        theta

DH_table = [0, -pi/2,     0.2,        0;...
            0, +pi/2,       0,        0;...
            0, -pi/2,     0.2,        0;...
            0, +pi/2,       0,        0;...
            0, -pi/2,     0.2,        0;...
            0, +pi/2,       0,        0;...
            0, -pi/2,     0.2,        0;...
            0, +pi/2,       0,        0;...
            0,     0,     0.2,        0];

Jtype_list(1:nj) = 'R';
Jtype_list(1) = 'P';
b_list = -[0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0, 0.2, 0]./2; % joint position w.r.t. DH frame

import casadi.*
q = SX.sym('q', nj, 1);  % symbolic variables initialization 

% symbolic forward kinematics with DH template transforms
[T0E, Tj, T0j] = DHFWkin(DH_table, q, Jtype_list);

% base offset, if any
Toffset0 = eye(4);
ToffsetE = eye(4);
Jac = DHJac(T0j, Jtype_list, 'baseOffset', Toffset0, 'eeOffset', ToffsetE);

% EE offset, if any
ToffsetE = eye(4);
T0E = T0E*ToffsetE;
Jac = twistPole(ToffsetE)*Jac;


% build functions from symbolic expressions
T0E_fun = Function('T0E', {q}, {T0E});
Tj_fun = Function('Tj', {q}, {Tj{:}});
Jac_fun = Function('Jac', {q}, {Jac});

%% definition of configuration errors
import casadi.*
Kp = 15;
Ko = 15;

% gain matrix
K = [SX.eye(3).*Kp, SX(3,3); SX(3,3), SX.eye(3).*Ko];

% orientation error with quaternion
[nT, thetaT] = rotToAxisAngle(R_expr);
[nE, thetaE] = rotToAxisAngle(T0E(1:3, 1:3));
Qt = axisAngleToQuat(thetaT, nT);
QE = axisAngleToQuat(thetaE, nE);

e_o = OerrorFromQuat(Qt, QE);

% orientation error with vecForm of skew sym matrix
% Rtilde = R_expr*T0E(1:3, 1:3).';
% e_o = vecForm(Rtilde);

% positional error
e_p = T_expr(1:3, 4)-T0E(1:3, 4);

% error vector
e = [e_p; e_o];

e_fun = Function('err', {t, q}, {e});

L = twistRefCLIK(T_expr, T0E);

I = SX.eye(3);
O = SX(3,3);
Linv = [I,      O;...
        O, inv(L)];
LT = [I, O;...
      O, L'];

%% PseudoInverse Jacobian and lower priority tasks

% avoid collision
radius_sphere = 0.1;
center_sphere = [0.1; -0.1; 0.3];
p_robot = Toffset0*T0j{6}(1:4, 4);
H = sum((p_robot(1:3) - center_sphere).^2); % ... + other joint collision checks and other collision spheres

% gradient of H
gradH = jacobian(H,q)';

% minimum norm solution weights
weights = 2*nj:-2:1;
weights = weights./sum(weights);

q0dot = 2.*gradH;
[Jr, b] = redundantInverseKin(Jac, weights, q0dot);

%% kinematics integration

% qdot_expr = Jr*Linv*(LT*twist_expr + K*e) + b; % ODE rhs with axis angle error orientation
qdot_expr = Jr*(twist_expr + K*e) + b; % ODE rhs with quaternion error

% integration with runge kutta
dt = 0.01;
t0 = ti;
tend = tf;

% set the initial configuration of the manipulator
q0 = 0+zeros(nj, 1);
q0(1) = 0;
q0(2) = -pi/4;
q0(4) = pi/4;
q0(6) = pi/4;
q0(8) = pi/4;
% q0 = [-2.72; -1.6843; -0.2843; 0.7370; -0.3826; 0.7763; 2.9671; 1.6225; 2.2424];
qsol = RK4(q, [], qdot_expr, dt, tend, q0, t0, [], 't_expr', t);

%% graphics

% compute numerical kinematics

T0E_num = full(T0E_fun(q0));
Tj_num = cell(nj, 1);
[Tj_num{:}] = Tj_fun(q0);
Tj_num = cellfun(@full, Tj_num, 'UniformOutput',false);


% manipulator graphics enclosed inside this function
transforms = init_plot_DH(Tj_num, b_list, Jtype_list, 'joint_len', 0.06, 'joint_r', 0.02, 'T0', Toffset0);

% ground graphics
plane_pt = [[1;1], [1;-1], [-1;-1], [-1;1]].*[0.5;0.7];
patch(gca, plane_pt(1,:), plane_pt(2,:), zeros(1, 4), 'facecolor', [.7 .7 .7]);

% writing board graphics
board_pt = Ttx(boardOrigin(1))*Tty(boardOrigin(2))*Ttz(boardOrigin(3))*TrotY(pi/2)*([[0.3;0.3;0;1], [0.3;-0.3;0;1], [-0.3;-0.3;0;1], [-0.3;0.3;0;1]].*[1;2;1;1]);
fill3(gca, board_pt(1,:), board_pt(2,:), board_pt(3, :), board_pt(3, :).*0, 'facecolor', '#454545');

% task line plot
line(writing_pts(1,:), writing_pts(2,:), writing_pts(3,:), 'color', 'w', 'linewidth', 1.2)

% plot fixed frame
plotFrame(eye(4), 'label', '0', 'scale', 0.7) % plot fixed frame

% plot task frame
task_transform = hgtransform(gca);
task_transform.Matrix = full(T_fun(0));
plotFrame(eye(4), 'label', 't', 'scale', 0.1, 'parent', task_transform); % plot task frame

% plot collision sphere
plotSphere(radius_sphere, center_sphere, 'parent', gca, 'FaceAlpha', 1, 'FaceColor', 'y');

% scene plot options
view(-45,45)
xlabel('x')
ylabel('y')
zlabel('z')

set(gca, 'zlim', [0 1]);
set(gca, 'xlim', [-0.7 0.7]);
set(gca, 'ylim', [-0.7 0.7]);
set(gca, 'box', 'on')

%% animation
tracked_line = line(nan(1, size(qsol, 2)),nan(1, size(qsol, 2)),nan(1, size(qsol, 2)), 'color', 'g', 'linewidth', 1.4);
for j = 1:size(qsol, 2)
    
   updatePlot_DH(Tj_fun, qsol(1:nj, j), Toffset0, ToffsetE, transforms)  % update manipulator joints graphics
   T0E_num = full(T0E_fun(qsol(1:nj, j)));           % compute E-E matrix
   task_transform.Matrix = full(T_fun(qsol(end, j)));% update task frame graphics

   % append numerical values for the actual E-E trajectory
   tracked_line.XData(j) = T0E_num(1,4);
   tracked_line.YData(j) = T0E_num(2,4);
   tracked_line.ZData(j) = T0E_num(3,4);
   drawnow
end

%% configuration errors plot

% numerical conf. errors calculation
e_num = full(e_fun(qsol(end, :), qsol(1:nj, :)));
t_num = qsol(end, :);



% compute determinant of numerical jacobian to check for singularities
det_jac = nan(1, size(qsol, 2));

for j = 1:size(qsol, 2)
    Jac_num = full(Jac_fun(qsol(1:nj, j)));
    det_jac(j) = full(det(Jac_num*Jac_num'));
    
end

plot_orientation_errors(e_num, t_num, 'lineColor', 'b');


figure('color', 'w');
plot(t_num, det_jac, 'linewidth', 1.4)
xlabel('t [s]')
ylabel('$det(JJ^T)$')
set(gca, 'Fontsize', 30);

plot_joint_trajectories(qsol, 'lineColor', 'b');