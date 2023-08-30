function J = bodyJac_ZYX(alpha)
%
% body jacobian of the ZYX Euler angles parametrization (maps omega_body to angles derivatives)
%

theta = alpha(2);
phi = alpha(3);

sphi = sin(phi);
cphi = cos(phi);
ctheta = cos(theta);
stheta = sin(theta);

J = [     -stheta,     0, 1;...
      ctheta*sphi,  cphi, 0;...
      cphi*ctheta, -sphi, 0];

end