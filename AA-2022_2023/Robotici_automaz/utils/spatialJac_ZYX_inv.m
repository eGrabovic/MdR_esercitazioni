function J = spatialJac_ZYX_inv(alpha)
%
% inverse jacobian of the Euler angle ZYX parametrization.
% map from twist to derivatives of parametrization: alpha_dot
%
psi = alpha(1);
theta = alpha(2);

spsi = sin(psi);
cpsi = cos(psi);
ctheta = cos(theta);
stheta = sin(theta);

J = [(cpsi*stheta)/ctheta, (spsi*stheta)/ctheta, 1;...
                    -spsi,                 cpsi, 0;...
              cpsi/ctheta,          spsi/ctheta, 0];

end