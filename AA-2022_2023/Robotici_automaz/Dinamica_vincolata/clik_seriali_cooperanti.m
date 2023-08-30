%% definizione task per c.r.
t_sym = casadi.MX.sym('t');
t0 = 0;
tf = 15;
p0 = [0;0;0.2];
pf = [0;0.15;0.55];
phi0 = 0;
phif = pi/2;

d = (t_sym-t0)./(tf-t0) ; % parametro temporale normalizzato
p_t = (pf-p0)*d + p0;
phi_t = (phif-phi0)*d + phi0;

qR_t = [p_t;0;0; phi_t];
qR_t_fun = casadi.Function('qR', {t_sym}, {qR_t});
T0R_task = T0R_fun(qR_t);
T0R_task_fun = casadi.Function('T0R_t', {t_sym}, {T0R_task});

%% animazione di ciò che vogliamo far fare al c.r.per scommentare selezionare codice e premere CTRL + T

% figure('color', 'w'); hold on; axis equal
% 
% % rigid body plot
% transformR = hgtransform('parent', gca);
% transformR.Matrix = full(T0R_task_fun(0));
% [pointsR, facesR] = createParallelepiped(t, [L/2, L/2], [0;1;0], 0);
% patch('faces', facesR, 'vertices', pointsR.',  'facecolor', 'green', 'Parent', transformR);
% plotFrame(eye(4), 'scale', 0.1, 'label', 'R', 'parent', transformR);
% 
% common_scene_options;
% 
% for i = 0:0.01:tf
%     transformR.Matrix = full(T0R_task_fun(i));
%     pause(0)
% end

%% calcolo twist desiderato
T0R_dot_task = reshape(jacobian(T0R_task, t_sym), 4, 4);
v_task = T0R_dot_task(1:3, 4);
R0R_task = T0R_task(1:3,1:3);
omega_task = vecForm(T0R_dot_task(1:3, 1:3)*R0R_task); % R^T*Rdot
twistR_task = [v_task; omega_task];

v_dot = jacobian(v_task, t_sym);
omega_dot = jacobian(omega_task, t_sym);
twistR_dot_task = [v_dot; omega_dot];

%% calcolo della cinematica differenziale inversa con CLIK per valutare le leggi orarie necessarie sui giunti
import casadi.*
Kp = 15;
Ko = 15;

% gain matrix
K = [eye(3).*Kp, zeros(3, 3); zeros(3, 3), eye(3).*Ko];
K = blkdiag(K, K);

% orientation error with quaternion   CTRL + R comment highlighted code ; CTRL + T un-comment highlighted code
[nT, thetaT] = rotToAxisAngle(R0R_task);
[nEA, thetaEA] = rotToAxisAngle(T0EA(1:3, 1:3));
[nEB, thetaEB] = rotToAxisAngle(T0EB(1:3, 1:3));
Qt = axisAngleToQuat(thetaT, nT);
QEA = axisAngleToQuat(thetaEA, nEA);
QEB = axisAngleToQuat(thetaEB, nEB);

e_oA = OerrorFromQuat(Qt, QEA);
e_oB = OerrorFromQuat(Qt, QEB);

% positional error
pA_t = T0R_task*[0;-L/2;0;1];
pB_t = T0R_task*[0; L/2;0;1];
e_pA = pA_t(1:3)-T0EA(1:3, 4);
e_pB = pB_t(1:3)-T0EB(1:3, 4);

% error vector
eA = [e_pA; e_oA];
eB = [e_pB; e_oB];
e = [eA; eB];

e_fun = Function('err', {t_sym, q}, {e});

weights = ones(njA+njB, 1);
q0dot = zeros(njA+njB, 1);
HJ = H*J;
[Jr, b] = redundantInverseKin(HJ, weights, q0dot);
Grasp_fun = casadi.Function('grasp', {qR}, {Grasp});

Htwist = H*((Grasp_fun(qR_t)'*twistR_task + K*e));

qdot = Jr*Htwist + b;

% integration with runge kutta
dt = 0.01;
t0 = 0;
tf = tf;

% parametri integrazione

q0A = zeros(njA, 1)+0.5;
q0B = zeros(njB, 1)-0.5;
q0R = [p0; 0; 0; phi0];

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
problem = struct('x', [qA; qB], 'p', qR, 'f', 0, 'g', err);

% costruzione problema di ottimizzazione. Vengono calcolate in automatico
% le derivate per valutare hessiano e gradiente del Lagrangiano. Viene
% stabilita un' interfaccia col solutore IPOPT (algoritmo interior-point, stato dell'arte)
S = casadi.nlpsol('S', 'ipopt', problem);

% chiamata al solutore con initial guess numeriche e bound sui vincoli
sol = S('x0', [q0A; q0B], 'p', q0R, 'ubg', 0, 'lbg', 0); % upper e lower bound 0 sui vincoli significa che sono equality

% estrazione della soluzione
q0 = full(sol.x);
q0A = q0(1:njA);
q0B = q0(njA+1:njA+njB);

q0 = [q0A; q0B];
q = [qA; qB];
qsol = RK4(q, [], qdot, dt, tf, q0, t0, [], 't_expr', t_sym);

%% interpolazione leggi giunti per calcolare velocità ed accelerazioni
import casadi.* 
q_des = MX(njA+njB, 1);
q_d_des = MX(njA+njB, 1);
q_dd_des = MX(njA+njB, 1);

% interpolazione con spline casadi
t_num = full(qsol(end,:));
t_num(end) = tf+1e-8;
for i = 1:njA+njB
    q_interp = casadi.interpolant('interp', 'bspline', {t_num}, qsol(i,:));
    q_des(i) = q_interp(t_sym);
    q_d_des(i) = jacobian(q_des(i), t_sym);
    q_dd_des(i) = jacobian(q_d_des(i), t_sym);
end

q_des_fun = casadi.Function('qd', {t_sym}, {q_des});
q_d_des_fun = casadi.Function('qd_d', {t_sym}, {q_d_des});
q_dd_des_fun = casadi.Function('qd_dd', {t_sym}, {q_dd_des});