function [T] = DH_transform(DH_par, q, Jtype)
%
% [a, alpha, d, theta]
%

% extracting DH parameters
a = DH_par(1);
alpha = DH_par(2);

if strcmpi(Jtype, 'P') % prismatic joint
    d = DH_par(3) + q;
    theta = DH_par(4);
elseif strcmpi(Jtype, 'R') % revolute joint
    d = DH_par(3);
    theta = DH_par(4)+ q;
end

% writing commong trig expressions
cth = cos(theta);
sth = sin(theta);
calf = cos(alpha);
salf = sin(alpha);

T = [cth, -calf*sth,  salf*sth,    a*cth;...
    sth,  calf*cth, -salf*cth,    a*sth;...
    0,      salf,      calf,        d;...
    0,         0,         0,        1];


end