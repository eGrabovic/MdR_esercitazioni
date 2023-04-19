function Jac = DHJac_CG(T0j, Jtype_list, cg_vec, body_index, options)
%
%
%

arguments
    T0j
    Jtype_list
    cg_vec = []; % center of mass (of i-th body) expressed in {S0} frame
    body_index = [];
    options.baseOffset = eye(4);
    options.eeOffset = eye(4);
end

cl = class(T0j{1});
Jac = zeros(6, length(T0j), cl);
n = body_index;

if strcmpi(Jtype_list(1), 'P')
    Jac(1:3, 1) = [0;0;1];
elseif strcmpi(Jtype_list(1), 'R')
    Jac(:, 1) = [cross([0;0;1], cg_vec); [0;0;1]];
end

if n == 1
    return
end

for j = 2:n
    
    T0j_1 = T0j{j-1};
    k = T0j_1(1:3, 3); % z axis of previous frame
    O_0j = T0j_1(1:3, 4); % origin of previous frame

    O_j_1C = cg_vec - O_0j;


    if strcmpi(Jtype_list(j), 'P')
        Jac(1:3, j) = k;
    elseif strcmpi(Jtype_list(j), 'R')
        Jac(:, j) = [cross(k, O_j_1C); k];
    end


end

Jac = twistTransf(options.baseOffset)*Jac;
if n == length(T0j)
    Jac = twistPole(options.eeOffset)*Jac;
end

end