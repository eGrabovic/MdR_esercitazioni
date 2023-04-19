function x_num = collocation_integrator(x_expr, u_expr, xdot_expr, dt, t_end, x0, t0, u_in, options)
%
% x_num = RK4(x_expr, u_expr, xdot_expr, dt, t_end, x0, t0, u_in, options)
%
% 4-th order Runge-kutta integrator applied to dynamic system in state
% space:
% x_expr    = state space variable vector expressed as casadi symbolic
% u_expr    = input variable vector expressed as casadi symbolic; if no
%             input set it to empty [] 
% xdot_expr = f(x_expr, u_expr, t_expr), the RHS of the dynamic equation 
%             expressed as casadi symbolic expression
%             
% dt        = time interval (numeric)
% t_end     = time horizon (numeric)
% x0        = initial state (numeric)
% t0        = initial time (numeric)
% u_in      = discrete system input over the time steps 
%             (numeric, casadi symbolic or empty)
% 
% options:
%
%  'Nsteps'   (default = 1) sets the number of integration steps per
%             integration interval
%
%  't_expr'   (default = empty variable) sets a casadi expression for an
%             explicit time dependency of the dynamic system which can
%             directly appear in the xdot expression (time variant system);
%             else it is considered time Invariant

arguments
    x_expr
    u_expr
    xdot_expr
    dt
    t_end
    x0
    t0
    u_in
    options.NstepsInterval = 1;
    options.t_expr = [];
    options.degree = 3;
end

% some error checks
if size(u_expr, 1) ~= size(u_in, 1)
    error('dimension mismatch between the casdadi expression for the system input ''u_expr'' and the discrete input array u_in');
end

% extract optional name-valued arguments
Nsteps = options.NstepsInterval;
t = options.t_expr;

% check if the system is TV
if ~isempty(t)
    x_expr = vertcat(x_expr, t);
    x0(end+1, :) = t0;            % concatenate time to the state space vector
    xdot_expr = vertcat(xdot_expr, 1);
end

% compute the total step number
steps = floor((t_end - t0)./dt);

% initialize the discrete numeric states
x_num = nan(length(x0), steps+1); % inizialization numeric matrix of the state evolution
x_num(:, 1) = x0;                 % initial state

if isempty(u_in)                  % if the array remains empty the system has no inupt
    u_in = zeros(1, steps);
end

% initialize time dependency, if system is time invariant it wont affect
% the results 

import casadi.* % import casadi package to use the built-in Function 
xdot_fun = Function('xdot_fun', {x_expr, u_expr}, {xdot_expr}); % convert to casadi function
solver = []; % initialize an empty solver (the solver is built at the first time step inside the Euler_BS call)
colloc = [];
% integrate over the time horizon
for i = 2:steps+1
    
%     x_num(:, i) = Direct_collocation_num_matlab(x_num(:, i-1), u_in(:, i-1), xdot_fun, dt, 'degree', options.degree); % carry the solver built once through the loop
    [x_num(:, i), solver, colloc] = Direct_collocation_num(x_num(:, i-1), u_in(:, i-1), xdot_fun, dt, 'degree', options.degree, 'Solver', solver, 'x_colloc_guess', colloc);
end

end