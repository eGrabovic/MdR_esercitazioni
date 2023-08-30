function [Jac, Jac_dq] = DHJac_CG_derivative(T0j, Jtype_list, cg_vec_body, body_index, options)
%
%
%

arguments
    T0j
    Jtype_list
    cg_vec_body = []; % center of mass (of i-th body) expressed in {S0} frame
    body_index = [];
    options.baseOffset = eye(4);
    options.eeOffset = eye(4);
    options.derivative = [];
end

cl = class(T0j{1});
Jac = zeros(6, length(T0j), cl);
n = body_index;
Jac_dq = cell(1, n);
T0j_dq = options.derivative;

cg_vec = T0j{body_index}(1:3,1:3)*cg_vec_body + T0j{body_index}(1:3,4); % centro di massa in {S0}

if strcmpi(Jtype_list(1), 'P')
    Jac(1:3, 1) = [0;0;1];
elseif strcmpi(Jtype_list(1), 'R')
    Jac(:, 1) = [cross([0;0;1], cg_vec); [0;0;1]];
end

if nargout > 1 && ~isempty(T0j_dq)% compute also derivatives
         for ff = 1:n
            Jac_dq{ff} = zeros(6, length(T0j), cl);
         end
    for kk = 1:n
        O_j_1C_dq = T0j_dq{body_index}{kk}(1:3, 1:3)*cg_vec_body + T0j_dq{body_index}{kk}(1:3, 4);
        if strcmpi(Jtype_list(1), 'R')
            Jac_dq{kk}(1:3, 1) = hat([0;0;1])*O_j_1C_dq;
        end
    end
end

if n == 1
    return
end

for j = 2:n
    
    T0j_1 = T0j{j-1};
    k = T0j_1(1:3, 3); % z axis of previous frame
    O0j = T0j_1(1:3, 4); % origin of previous frame

    Oj_1C = cg_vec - O0j;


    if strcmpi(Jtype_list(j), 'P')
        Jac(1:3, j) = k;
    elseif strcmpi(Jtype_list(j), 'R')
        Jac(:, j) = [cross(k, Oj_1C); k];
    end
    
    if nargout > 1 && ~isempty(T0j_dq)% compute also derivatives
        for kk = 1 : n
            if kk >= j % prevuois joint does not depend on current joint variable
                dk_dq = [0;0;0];
                Oj_1_dq = [0;0;0];
            else
                dk_dq = T0j_dq{j-1}{kk}(1:3, 3);
                Oj_1_dq = T0j_dq{j-1}{kk}(1:3, 4);
            end
            O_j_1C_dq = T0j_dq{n}{kk}(1:3, 1:3)*cg_vec_body + T0j_dq{n}{kk}(1:3, 4) - Oj_1_dq;
            if strcmpi(Jtype_list(j), 'P')
                Jac_dq{kk}(1:3, j) = dk_dq;
            elseif strcmpi(Jtype_list(j), 'R')
                Jac_dq{kk}(:, j) = [hat(dk_dq)*Oj_1C+hat(k)*O_j_1C_dq; dk_dq];
            end
            
        end
    end
end

Jac = twistTransf(options.baseOffset)*Jac;
if n == length(T0j)
    Jac = twistPole(options.eeOffset)*Jac;
end
if nargout > 1 && ~isempty(options.derivative)% compute also derivatives
    for kk = 1 : n
        Jac_dq{kk} = twistTransf(options.baseOffset)*Jac_dq{kk};
        if n == length(T0j)
            Jac_dq{kk} = twistPole(options.eeOffset)*Jac_dq{kk};
        end
    end
end

end


