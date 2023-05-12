function Jac = DHJac(T0j, Jtype_list, options)
%
%
%

arguments
    T0j
    Jtype_list
    options.joint_index = [];
    options.baseOffset = eye(4);
    options.eeOffset = eye(4);
end

n = length(T0j);
if ~isempty(options.joint_index)
    n = options.joint_index;
end

O_0n = T0j{n}(1:3, 4);
cl = class(T0j{1});
Jac = zeros(6, length(T0j), cl);


if strcmpi(Jtype_list(1), 'P')
    Jac(3, 1) = 1;
elseif strcmpi(Jtype_list(1), 'R')
    Jac(:, 1) = [cross([0;0;1], O_0n); [0;0;1]];
end

for j = 2:n

    T0j_1 = T0j{j-1};
    k = T0j_1(1:3, 3);
    O_0j = T0j_1(1:3, 4);

    O_jn = O_0n - O_0j;


    if strcmpi(Jtype_list(j), 'P')
        Jac(1:3, j) = k;
    elseif strcmpi(Jtype_list(j), 'R')
        Jac(:, j) = [cross(k, O_jn); k];
    end


end

Jac = twistTransf(options.baseOffset)*Jac;
if n == length(T0j)
    Jac = twistPole(options.eeOffset)*Jac;
end
end