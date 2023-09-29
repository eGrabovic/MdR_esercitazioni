function J = bodyJac_ZYX_inv(alpha)
%
% inverse of body jacobian of the ZYX Euler angles parametrization (maps angles derivatives to body angular vel)
%

theta = alpha(2);
phi = alpha(3);

sphi = sin(phi);
cphi = cos(phi);
ctheta = cos(theta);
stheta = sin(theta);

J = [0,          sphi/ctheta,          cphi/ctheta;...
     0,                 cphi,                -sphi;...
     1, (sphi*stheta)/ctheta, (cphi*stheta)/ctheta];

end