%% Esercitazione CLICK 9R
% 
% Analisi cinematica differenziale inversa con inseguimento di traiettorie
% con algoritmo CLIK (Closed-Loop Inverse Kinematics)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex'); 
set(0,'defaultTextInterpreter','latex');

clc; clear all;

addpath(genpath('../graphics'))
addpath(genpath('../utils'))
addpath('Data')

%% definizione task

% import coordinates for the writing task
load('coordinateX_scritta.mat'); % variabile X
load('coordinateY_scritta.mat'); % variabile Y

m = 100;  % number of DFT coefficients

% discrete Fourier transform
pk = DFT_vec([X;Y], m);

% inverse Fourier transform (the actual analytic expression of the trajectory)
p = @(t) IFST(t, pk, m);

Nd = 2000; % number of sampling points

alpha = linspace(0,2*pi,Nd); % discretized time period

norm_pk = vecnorm(pk, 2, 1);
[norm_pk, idorder] = sort(norm_pk, 'descend');
%% animation of the curve
% Nd = 2000; % number of sampling points
% 
% alpha = linspace(0,2*pi,Nd); % discretized time period
% 
% norm_pk = vecnorm(pk, 2, 1);
% [norm_pk, idorder] = sort(norm_pk, 'descend');
% 
% figure; hold on; axis equal
% pl = plot(nan(1,Nd), nan(1,Nd));
% 
% pl_k = gobjects(length(pk), 1);
% for i = 1:length(pk)
%    pl_k(i) =  plot(nan(1,2), nan(1,2), 'color', 'k', 'linewidth', 1.4);
% end
% 
% for i = 1:Nd
%     
%    [pt, pt_k, pt_kp] = IFST(alpha(i), pk, m, idorder);
%    pl.XData(i)  = pt(1);
%    pl.YData(i)  = pt(2);
%    
%    for k = 1:length(pt_k)
%        pl_k(k).XData = [pt_kp(1,k), pt_k(1,k)];
%        pl_k(k).YData = [pt_kp(2,k), pt_k(2,k)];
%    end
%    drawnow
% end

%% Task trajectory calculation 
import casadi.*
boardOrigin = [0.35;0;0.5];         % center of the writing board
t = SX.sym('t', 1, 1);                    % symbolic variable for time
ti = 0;                             % initial time
tf = 30;                            % final time
alpha = (t-ti)./(tf-ti).*2*pi;      % scaling of the Fourier transofrm period
pt_expr = IFST(alpha, pk, m, idorder).*0.0015; % scaling factor to resize the writing
pt_expr = Ttransl(boardOrigin)*TrotY(-pi/2)*TrotZ(3*pi/2)*vertcat(pt_expr, 0, 1); % re-orienting the writing on the board
pt_fun = Function('pt', {t}, {pt_expr});
writing_pts = nan(4, Nd);
t_num = linspace(ti, tf, Nd);
for i = 1:Nd
    writing_pts(:, i) = full(pt_fun(t_num(i)));
end

%% calcolo del twist task
import casadi.*
% velocità lineare 
v_expr = jacobian(pt_expr, t); % der_pd/der_t

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

t_num = linspace(ti, tf, Nd);
twist_num = full(twist_fun(t_num));


figure('color', 'w')
ax = subplot(2,3,1);
plot(ax, t_num, twist_num(1,:))
xlabel('$t$ (s)');
ylabel('$v_x$ (m/s)');
ax = subplot(2,3,2);
plot(ax, t_num, twist_num(2,:))
xlabel('$t$ (s)');
ylabel('$v_y$ (m/s)');
ax = subplot(2,3,3);
plot(ax, t_num, twist_num(3,:))
xlabel('$t$ (s)');
ylabel('$v_z$ (m/s)');
ax = subplot(2,3,4);
plot(ax, t_num, twist_num(4,:))
xlabel('$t$ (s)');
ylabel('$\omega_x$ (m/s)');
ax = subplot(2,3,5);
plot(ax, t_num, twist_num(5,:))
xlabel('$t$ (s)');
ylabel('$\omega_y$ (m/s)');
ax = subplot(2,3,6);
plot(ax, t_num, twist_num(6,:))
xlabel('$t$ (s)');
ylabel('$\omega_z$ (m/s)');
set(findall(gcf,'-property','FontSize'),'FontSize',22)
sgtitle('Reference signals')

%% redundant robot parametrization

% robot 9R
nj = 9; % number of joints
d = 0.10;
z = repmat(d, 1, nj);
z(1) = 0;

% creation of offset transformations
g0 = cell(1,nj);
X = zeros(6, nj);
Jtype(1:nj) = 'R';
for ii = 1:nj
    
   g0{ii} = Ttz(z(ii));
   if mod(ii,2)
       X(:, ii) = [0;0;0;0;0;1]; % all revolute joints
   else
       X(:, ii) = [0;0;0;1;0;0]; % all revolute joints
   end
    
end
% Jtype(1) = 'P';
% X(:, 1) = [0;0;1;0;0;0];
g0{ii+1} = eye(4); % placeholder transform for graphical purposes (E-E drawing)

import casadi.*
q = casadi.SX.sym('q', nj, 1);  % symbolic variables initialization 
Toffset0 = eye(4);%Ttz(1)*TrotX(pi);
ToffsetE = Ttz(0.08);
[G0E, Grel] = FWkin_localPOE(g0, X, q);
G0E = Toffset0*G0E*ToffsetE;
Grel{1} = Toffset0*Grel{1};
Jac_expr = spatialJac_localPOE(g0, X, q, 'EEoffset', ToffsetE, 'Toffset0', Toffset0);

Jac_expr = twistPole(-G0E(1:3, 4))*Jac_expr;

% build functions from symbolic expressions
T0E_fun = Function('T0E', {q}, {G0E});
Tj_fun = Function('Tj', {q}, {Grel{:}});
Jac_fun = Function('Jac', {q}, {Jac_expr});

%% definition of configuration errors
import casadi.*
Kp = 15;
Ko = 15;

% gain matrix
K = [SX.eye(3).*Kp, SX(3,3); SX(3,3), SX.eye(3).*Ko];

% orientation error with quaternion   CTRL + R comment highlighted code ; CTRL + T un-comment highlighted code
[nT, thetaT] = rotToAxisAngle(R_expr);
[nE, thetaE] = rotToAxisAngle(G0E(1:3, 1:3));
Qt = axisAngleToQuat(thetaT, nT);
QE = axisAngleToQuat(thetaE, nE);

e_o = OerrorFromQuat(Qt, QE);

% orientation error with vecForm of skew sym matrix
% Rtilde = R_expr*G0E(1:3, 1:3).';
% e_o = vecForm(Rtilde - Rtilde')./2;

% positional error
e_p = T_expr(1:3, 4)-G0E(1:3, 4);

% error vector
e = [e_p; e_o];

e_fun = Function('err', {t, q}, {e});

L = twistRefCLIK(T_expr, G0E);

I = SX.eye(3);
O = SX(3,3);
Linv = [I,      O;...
        O, inv(L)];
LT = [I, O;...
      O, L'];
  
L_fun = Function('L', {t, q}, {L});

%% PseudoInverse Jacobian and lower priority tasks

% avoid singularities
H = sqrt(det(Jac_expr*Jac_expr.'));

% stay close to a specified configuration
% q0 = 0.05+zeros(nj, 1);
% q0(1) = 0;
% q0(2) = -pi/4;
% q0(4) = pi/4;
% q0(6) = pi/4;
% q0(8) = pi/4;
% qmin = -pi + q0;
% qmax = pi + q0;
% qbar = (qmax+qmin)./2;
% H = -1./(2*nj)*sum((q - qbar)./qmax-qmin);

% gradient of H
gradH = jacobian(H,q)';

% minimum norm solution weights
weights = 2*nj:-2:1;
weights(end) = 0.1;
weights = weights./sum(weights);
q0dot = 2'.*gradH;
[Jr, b] = redundantInverseKin(Jac_expr, weights, q0dot);

%% kinematics integration

% qdot_expr = Jr*Linv*(LT*twist_expr + K*e) + b; % ODE rhs with axis angle error orientation
qdot_expr = Jr*(twist_expr + K*e) + b; % ODE rhs with quaternion error

% integration with runge kutta
dt = 0.01;
t0 = ti;
tend = tf;

% set the initial configuration of the manipulator
q0 = pi/6+zeros(nj, 1);
q0(2) = -pi/4;

qsol = RK4(q, [], qdot_expr, dt, tend, q0, t0, [], 't_expr', t);

%% graphics

% manipulator graphics enclosed inside this function
[ax, transforms] = plotRobot_localPOE(Tj_fun, g0, X, q0, Jtype);


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
for j = 1:1:size(qsol, 2)
    
   updatePlot_localPOE(Tj_fun, qsol(1:nj, j), transforms);
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