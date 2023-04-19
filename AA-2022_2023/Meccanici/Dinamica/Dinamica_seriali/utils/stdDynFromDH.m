function [B, C, G, B_dq] = stdDynFromDH(DH_table, Jtype_list, q, qdot, cg_list, mass_list, inertia_list, options)
%
% input: - DH_table: DH parameters (each row = [a, alpha, d, theta])
%        - cg_list: center of mass list of each link (from origin of S_{i-1}) 
%        - inertia_list: generalized inertia of each link
%
arguments
    DH_table
    Jtype_list
    q
    qdot
    cg_list
    mass_list
    inertia_list
    options.K = 0;
    options.q_eq = 0;
    options.gravity = [0;0;-9.81];
    options.baseOffset = eye(4);
    options.eeOffset = eye(4);
end
% extract data type we are working with (double or casadi symbolic)
% and initialization of data
import casadi.*
cl = class(q);
N = size(cg_list, 2);

B = zeros(N, N, cl);
C = zeros(N, N, cl);
G = zeros(N, 1, cl);

% forward kinematics map
[~, ~, T0j] = DHFWkin(DH_table, q, Jtype_list);
Toffset0 = options.baseOffset;

% compute cg jacobian and inertia matrix
for j = 1:N
    T0jOff = T0j{j};
    cg_vec = T0jOff*[cg_list(:, j); 1]; % posizione centro di massa in {S0}
    Jac = DHJac_CG(T0j, Jtype_list, cg_vec(1:3), j, 'baseOffset', Toffset0, 'eeOffset', options.eeOffset);
    T0jOff = Toffset0*T0jOff;
    inertia = T0jOff(1:3,1:3)*inertia_list{j}*(T0jOff(1:3,1:3)');
    Jv = Jac(1:3, :);
    Jw = Jac(4:6, :);
    B = B + mass_list(j).*Jv'*Jv + Jw'*inertia*Jw; % sommatoria di B_i
    G = G - mass_list(j).*(Jv'*options.gravity);                                              % sommatoria di G_i
end


B_dq = cell(1, N);

B_dq_sym = jacobian(B, q);
for i = 1:N
    B_dq{i} = reshape(B_dq_sym(:, i), N, N);
end

G = G - options.K*(q-options.q_eq); % add elastic forces contribution

% calcolo matrice C coi simboli di Christoffel del primo tipo
for i = 1:N
    for j = 1:N
        Bij_dqk = zeros(1, N, cl);
        for k = 1:N
            Bij_dqk(k) = B_dq{k}(i, j);
        end
        Bik_dqj = B_dq{j}(i, :);
        Bjk_dqi = B_dq{i}(j, :);
        C(i, j) = C(i, j) + 0.5.*(Bij_dqk + Bik_dqj - Bjk_dqi)*qdot;
    end
end

end