%% DH table
nj = 6;

           %a     %alpha     %d    %theta
DH_table = [0, -pi/2,  0.4, 0;...
            0, pi/2,  0.4, 0;...
            0, 0,  0, 0;...
            0, pi/2,  0.4, -pi/2;...
            0, pi/2,  0, +pi/2;...
            0, 0,  0.4, 0];

Jtype_list(1:nj) = 'R'; Jtype_list(3) = 'P';

q = casadi.SX.sym('q', 6, 1);

[T0E, Tj, T0j] = DHFWkin(DH_table, q, Jtype_list);

% build functions from symbolic expressions
T0E_fun = casadi.Function('T0E', {q}, {T0E});
Tj_fun = casadi.Function('Tj', {q}, {Tj{:}});
Toffset0 = eye(4);ToffsetE = eye(4);

b_list = zeros(1, 6);
b_list(4) = -0.2;
b_list(6) = -0.2;

%% graphics
q0 = zeros(6,1);
T0E_num = full(T0E_fun(q0));
Tj_num = cell(nj, 1);
[Tj_num{:}] = Tj_fun(q0);
Tj_num = cellfun(@full, Tj_num, 'UniformOutput', false);

transforms = init_plot_DH(Tj_num, b_list, Jtype_list, 'joint_len', 0.06, 'joint_r', 0.02, 'T0', Toffset0);
q0(3) = 0.05;
updatePlot_DH(Tj_fun, q0, Toffset0, ToffsetE, transforms)