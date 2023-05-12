function [T0E, Tj, T0j] = DHFWkin(DH_table, q, Jtype_list)
%
% [T0E, Tj, T0j] = DHFWkin(DH_table, q, Jtype_list)
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
% OUTPUT:
% - T0E: from ground to E-E transformation
% - Tj: from j-1 to j transformation
% - T0j: from ground to j transformation
%
arguments
    DH_table
    q
    Jtype_list
end

% 
Tj = cell(1, length(q));
T0j = cell(1, length(q));
% first joint
Tj{1} = DH_transform(DH_table(1, :), q(1), Jtype_list(1));
T0E = Tj{1};
T0j{1} = Tj{1};

% from second to E-E
for j=2:length(q)

    Tj{j} = DH_transform(DH_table(j, :), q(j), Jtype_list(j));
    T0E = T0E*Tj{j};
    T0j{j} = T0E;

end

end

