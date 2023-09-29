function [T, dT_dq] = DH_transform(DH_par, q, Jtype)
%
% [a, alpha, d, theta]
%

% extracting DH parameters
a = DH_par(1);
alpha = DH_par(2);

if strcmpi(Jtype, 'P') % prismatic joint
    d = DH_par(3) + q;
    d_dq = 1;
    theta = DH_par(4);
    dtheta_dq = 0;
elseif strcmpi(Jtype, 'R') % revolute joint
    d = DH_par(3);
    d_dq = 0;
    theta = DH_par(4)+ q;
    dtheta_dq = 1;
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

if nargout > 1
    dT_dq = [-sth,     -cth*calf,     cth*salf,     -a*sth;...
              cth,     -sth*calf,     sth*salf,      a*cth;...
                0,             0,            0,          0;...
                0,             0,            0,          0].*dtheta_dq...
                +...
                [0, 0, 0,    0;...
                 0, 0, 0,    0;...
                 0, 0, 0, d_dq;
                 0, 0, 0,    0];
end

end