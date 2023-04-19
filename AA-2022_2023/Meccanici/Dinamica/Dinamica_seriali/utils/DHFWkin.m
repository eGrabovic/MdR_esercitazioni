function [T0E, Tj, T0j, Tj_dq] = DHFWkin(DH_table, q, Jtype_list)
%
% [T0E, Tj] = DHFWkin(DH_table, q, Jtype_list)
%
% forward kinematics calculation with DH parameterization
% 
% DH_table = DH parameters (each row = [a, alpha, d, theta])
%
% q = joint values (elongation/rotation)
%
% Jtype_list = array of strings which defines the joiny type:
%               'P' for prismatic
%               'R' for revolute
%
arguments
    DH_table
    q
    Jtype_list
end

if nargout > 3 % output with differentiation (extra computations)
    Tj = cell(1, length(q));
    T0j = cell(1, length(q));
    Tj_dq = cell(1, length(q));
    [Tj{1}, Tj_dq{1}] = DH_transform(DH_table(1, :), q(1), Jtype_list(1));
    T0E = Tj{1};
    T0j{1} = Tj{1};
    
    for j=2:length(q)
        
        [Tj{j}, Tj_dq{j}] = DH_transform(DH_table(j, :), q(j), Jtype_list(j));
        T0E = T0E*Tj{j};
        T0j{j} = T0E;
        
    end
    
    return
    
end

% output without differentiation
Tj = cell(1, length(q));
T0j = cell(1, length(q));
Tj{1} = DH_transform(DH_table(1, :), q(1), Jtype_list(1));
T0E = Tj{1};
T0j{1} = Tj{1};

for j=2:length(q)

    Tj{j} = DH_transform(DH_table(j, :), q(j), Jtype_list(j));
    T0E = T0E*Tj{j};
    T0j{j} = T0E;

end

end

